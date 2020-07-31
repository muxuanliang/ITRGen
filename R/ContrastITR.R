ContrastITR <- function(data, is.weight = TRUE, x.test = NULL, nfold = 10){
  library(dplyr)
  size <- dim(data$predictor)[1]
  if(is.weight){
    density.ratio <- densratio::densratio(data$predictor, x.test)
    sample.weight <- pmin(1/density.ratio$compute_density_ratio(data$predictor), rep(10^4, times=size))
  }
  fit <- NULL
  fit_tree <- grf::causal_forest(X=data$predictor, Y=data$outcome, W=data$treatment,tune.parameters = 'all')
  contrast <- predict(fit_tree)$predictions * sample.weight

  fit <- DRITR(data=list(X=data$predictor, C=contrast), CONST = 1, maxit.mm = 100)
  fit
}

library(dplyr)
library(lbfgs)
DCA <- function(    # CONST = 1
  data, par.init = NULL, cnst.l2 = Inf, pen.l2 = 1e-8,
  method.mm = "proj_gd", tol.gd.mm = tol.gd, tol.obj.mm = tol.obj, tol.par.mm = tol.par, maxit.mm = maxit,
  tol.gd = 1e-5, tol.obj = 1e-5, tol.par = 1e-5, maxit = 200, verbose = F)
{
  n <- nrow(data$X)
  sv.X1 <- with(data, svd(cbind(1,X))$d[1])
  Lips.upper0 <- 2/n * svd(with(data, t(cbind(1,X)) %*% diag(abs(C)) %*% cbind(1,X)))$d[1]^2 + pen.l2

  par <- c(par.init[1], proj.l2(par.init[-1], cnst.l2))
  for(iter in 1:maxit) {
    # if(!verbose) {
    #   print(paste0("The ", iter, "-th DCA"))
    # }

    slope.l2 <- sqrt(sum(par[-1]^2))
    pen <- pen.l2/2 * slope.l2^2
    pen.gd <- pen.l2 * c(0, par[-1])
    obj <- OBJ(par, data) + pen
    obj.gd <- OBJ.gd(par, data) + pen.gd
    obj.gd.l2 <- sqrt(sum(obj.gd^2))

    f <- with(data, par[1] + as.vector(X %*% par[-1]))
    data$S <- with(data, pmax(C, 0) * surgN.gd(f) + pmax(-C, 0) * surgP.gd(f))

    if(verbose) {
      print(paste0("The ", iter, "-th iter: l2-norm of the slope=", slope.l2,
                   "; obj=", obj,
                   "; l2-norm of obj.gd=", obj.gd.l2))
    }

    if(obj.gd.l2 < tol.gd * max(slope.l2,1)) {
      message("Majorant update: gradient converges!"); break
    }
    # if(iter >= 2) {
    #   if(abs(obj - obj_old) < tol.obj * abs(obj_old)) {
    #     message("Majorant update: objective converges!"); break
    #   }
    # }

    par_old <- par
    obj_old <- obj
    if(method.mm == "proj_gd") {
      for(iter.mm in 1:maxit.mm) {
        par0 <- par
        gd0 <- L_ERM.gd(par0, pen.l2, data)
        if(sum(gd0^2) < tol.gd.mm^2 * max(sum(par0^2),1)) break
        L0 <- L_ERM(par0, pen.l2, data)
        f0 <- with(data, par0[1] + as.vector(X %*% par0[-1]))
        W0 <- with(data, abs(C) * (C > 0 & f0 >= 0 & f0 <= 1 | C < 0 & f0 >= -1 & f0 <= 0))
        Lips.lower <- 2/n * sv.X1^2 * mean(W0) + pen.l2
        Lips.upper <- min(Lips.upper0, 2/n * sv.X1^2 * max(W0) + pen.l2)
        ### line-search
        Lips.ls <- Lips.lower
        repeat {
          par1 <- par0 - gd0/Lips.ls  # one-step gradient descent
          par <- c(par1[1], proj.l2(par1[-1], cnst.l2))
          if(L_ERM(par, pen.l2, data) <= L0 + sum(gd0 * (par - par0)) + Lips.ls/2 * sum((par - par0)^2)) {
            break  # Armijo condition: successful majorization
          } else if(Lips.ls >= Lips.upper) {
            warning(paste0("The ", iter, "-th majorant minimization: gradient descent failed to majorize!"))
            break
          }
          Lips.ls <- min(2 * Lips.ls, Lips.upper)
        }
        if(sum((par - par0)^2) <= tol.par.mm^2 * max(sum(par0^2),1)) break
      }
      if(iter.mm == maxit.mm) {
        warning(paste0("The ", iter, "-th majorant minimization: maximal iteration is attained!"))
      } else if(verbose) {
        print(paste0("Equivalent l2-penalty parameter=",
                     pen.l2_equiv <- max(sqrt(sum(L_ERM.gd(par,0,data)[-1]^2) / sum(par[-1]^2)), pen.l2)))
      }
    } else if(method.mm == "lbfgs") {
      pen.l2_equiv <- pen.l2
      for(ftol in 1/2 * 10^seq(0, -4, -1)) {
        optim <- lbfgs(
          L_ERM, L_ERM.gd, par, pen.l2 = pen.l2, data = data, epsilon = tol.gd.mm, delta = tol.obj.mm,
          linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING", ftol = ftol, past = 20,
          max_iterations = maxit.mm, invisible = 1)
        if(optim$convergence != -998) break
      }
      if(optim$convergence == -998) {
        warning(paste0("The ", iter, "-th majorant minimization: maximal iteration is attained!"))
      } else if(ftol != 1/2) {
        warning(paste0("The ", iter, "-th majorant minimization: line-search accuracy parameter ", ftol))
      }
      if(optim$convergence == -997) {
        warning(paste0("The ", iter, "-th majorant minimization: L-BFGS doesn't converge!"))
      }
      par <- optim$par
    } else {
      stop('method.mm = "proj_gd" or "lbfgs"')
    }
    if(sum((par - par_old)^2) < tol.par^2 * max(sum(par_old^2),1)) {
      message("Majorant update: parameter converges!"); break
    }
  }
  if(iter == maxit) {
    warning("Majorant update: maximal iteration attained!")
  }
  return(list(par = par, eta = NA, lambda = NA, pen.l2 = ifelse(exists("pen.l2_equiv"), pen.l2_equiv, NA)))
}
PrDCA <- function(  # CONST > 1, k = Inf
  data, par.init = NULL, cnst.l2 = Inf, pen.l2 = 1e-8, CONST = 1,
  method.mm = "proj_gd", tol.gd.mm = tol.gd, tol.obj.mm = tol.obj, tol.par.mm = tol.par, maxit.mm = maxit,
  tol.gd = 1e-5, tol.obj = 1e-5, tol.par = 1e-5, maxit = 200, tol.CVaR = 1e-8, verbose = F)
{
  n <- nrow(data$X)
  ALPHA <- 1 - 1/CONST
  sv.X1 <- with(data, svd(cbind(1,X))$d[1])
  z <- with(data, c(C, -C))
  z.lowsum <- numeric()  # z.lowsum[i] = sum(pmax(z[i] - z, 0))
  id.sort <- sort(z, index.return = T)$ix
  z.lowsum[id.sort] <- 1:length(z) * z[id.sort] - cumsum(z[id.sort])
  data$ZP <- CONST/2 * z.lowsum[1:n]
  data$ZN <- CONST/2 * z.lowsum[(n+1):(2*n)]

  par <- c(par.init[1], proj.l2(par.init[-1], cnst.l2))
  for(iter in 1:maxit) {
    # if(!verbose) {
    #   print(paste0("The ", iter, "-th enhanced PrDCA"))
    # }

    slope.l2 <- sqrt(sum(par[-1]^2))
    pen <- pen.l2/2 * slope.l2^2
    pen.gd <- pen.l2 * c(0, par[-1])
    obj <- OBJ(par, data, Inf, CONST) + pen
    obj.gd <- OBJ.gd(par, data, Inf, CONST) + pen.gd
    obj.gd.l2 <- sqrt(sum(obj.gd^2))

    f <- with(data, par[1] + as.vector(X %*% par[-1]))
    prob <- c(surg(f), surg(-f))
    eta.tol <- VaR(z, prob, ALPHA, tol.CVaR)
    CVaR.tol <- CVaR(eta.tol, z, prob, ALPHA)
    eta <- ifelse(length(eta.tol) == 1, eta.tol, sample(eta.tol, 1))
    eps <- CVaR.tol[eta.tol == eta] - min(CVaR.tol)
    ZP_eta <- with(data, CONST/2 * pmax(+C - eta, 0))
    ZN_eta <- with(data, CONST/2 * pmax(-C - eta, 0))
    data$S <- with(data, ZP * surgP.gd(f) - ZN * surgP.gd(-f) -
                     (ZP_eta * surg.gd(+f) - ZN_eta * surg.gd(-f)))

    if(verbose) {
      print(paste0("The ", iter, "-th iter: eta=", eta,
                   "; l2-norm of the slope=", slope.l2,
                   "; obj=", obj,
                   "; l2-norm of obj.gd=", obj.gd.l2))
    }

    if(obj.gd.l2 < tol.gd * max(slope.l2,1)) {
      message("Majorant update: gradient converges!"); break
    }
    # if(iter >= 2) {
    #   if(abs(obj - obj_old) < tol.obj * abs(obj_old)) {
    #     message("Majorant update: objective converges!"); break
    #   }
    # }

    par_old <- par
    obj_old <- obj
    L_old <- L(par_old, pen.l2, data)
    if(method.mm == "proj_gd") {
      for(iter.mm in 1:maxit.mm) {
        par0 <- par
        gd0 <- L.gd(par0, pen.l2, data)
        if(sum(gd0^2) < tol.gd.mm^2 * max(sum(par0^2),1)) break
        L0 <- L(par0, pen.l2, data)
        f0 <- with(data, par0[1] + as.vector(X %*% par0[-1]))
        W0 <- with(data, ZP * (f0 >= 0 & f0 <= 1) + ZN * (f0 >= -1 & f0 <= 0))
        Lips.lower <- 2/n * sv.X1^2 * mean(W0) + pen.l2
        Lips.upper <- 2/n * sv.X1^2 * max(W0) + pen.l2
        ### line-search
        Lips.ls <- Lips.lower
        repeat {
          par1 <- par0 - gd0/Lips.ls  # one-step gradient descent
          par <- c(par1[1], proj.l2(par1[-1], cnst.l2))
          if(L(par, pen.l2, data) <= L0 + sum(gd0 * (par - par0)) + Lips.ls/2 * sum((par - par0)^2)) {
            break  # Armijo condition: successful majorization
          } else if(Lips.ls >= Lips.upper) {
            warning(paste0("The ", iter, "-th majorant minimization: gradient descent failed to majorize!"))
            break
          }
          Lips.ls <- min(2 * Lips.ls, Lips.upper)
        }
        if(sum((par - par0)^2) <= tol.par.mm^2 * max(sum(par0^2),1)) break
      }
      if(iter.mm == maxit.mm) {
        warning(paste0("The ", iter, "-th majorant minimization: maximal iteration is attained!"))
      } else if(verbose) {
        print(paste0("Equivalent l2-penalty parameter=",
                     pen.l2_equiv <- max(sqrt(sum(L.gd(par,0,data)[-1]^2) / sum(par[-1]^2)), pen.l2)))
      }
    } else if(method.mm == "lbfgs") {
      pen.l2_equiv <- pen.l2
      for(ftol in 1/2 * 10^seq(0, -4, -1)) {
        optim <- lbfgs(
          L, L.gd, par, data = data, epsilon = tol.gd.mm, delta = tol.obj.mm,
          linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING", ftol = ftol, past = 20,
          max_iterations = maxit.mm, invisible = 1)
        if(optim$convergence != -998) break
      }
      if(optim$convergence == -998) {
        warning(paste0("The ", iter, "-th majorant minimization: L-BFGS maximum number of line-searches attained!"))
      } else if(ftol != 1/2) {
        warning(paste0("The ", iter, "-th majorant minimization: line-search accuracy parameter ", ftol))
      }
      if(optim$convergence == -997) {
        warning(paste0("The ", iter, "-th majorant minimization: L-BFGS doesn't converge!"))
      }
      par <- optim$par
    } else {
      stop('method.mm = "proj_gd" or "lbfgs"')
    }
    if(L_old - L(par, pen.l2, data) < eps) {
      par <- par_old
    } else if(sum((par - par_old)^2) < tol.par^2 * max(sum(par_old^2),1)) {
      message("Majorant update: parameter converges!"); break
    }
  }
  if(iter == maxit) {
    warning("Majorant update: maximal iteration attained!")
  }
  return(list(par = par, eta = eta, lambda = NA, pen.l2 = ifelse(exists("pen.l2_equiv"), pen.l2_equiv, NA)))
}
BSUM <- function(   # CONST > 1, 1 < k <= Inf
  data, par.init = NULL, cnst.l2 = Inf, pen.l2 = 1e-8, k = 2, CONST = 1,
  method.mm = "proj_gd", tol.gd.mm = tol.gd, tol.obj.mm = tol.obj, tol.par.mm = tol.par, maxit.mm = maxit,
  tol.gd = 1e-5, tol.obj = 1e-5, tol.par = 1e-5, maxit = 200, tol.eta = 1e-5, tol.CVaR = 1e-8, verbose = F)
{
  n <- nrow(data$X)
  ALPHA <- 1 - 1/CONST
  k_star <- ifelse(k == Inf, 1, k/(k-1))
  sv.X1 <- with(data, svd(cbind(1,X))$d[1])
  z <- with(data, c(C, -C))

  par <- c(par.init[1], proj.l2(par.init[-1], cnst.l2))
  for(iter in 1:maxit) {
    # if(!verbose) {
    #   print(paste0("The ", iter, "-th BSUM"))
    # }

    slope.l2 <- sqrt(sum(par[-1]^2))
    pen <- pen.l2/2 * slope.l2^2
    pen.gd <- pen.l2 * c(0, par[-1])
    obj <- OBJ(par, data, k, CONST) + pen
    obj.gd <- OBJ.gd(par, data, k, CONST) + pen.gd
    obj.gd.l2 <- sqrt(sum(obj.gd^2))

    f <- with(data, par[1] + as.vector(X %*% par[-1]))
    prob <- c(surg(f), surg(-f))
    if(k == Inf) {
      eta.tol <- VaR(z, prob, ALPHA, tol.CVaR)
      CVaR.tol <- CVaR(eta.tol, z, prob, ALPHA)
      eta <- ifelse(length(eta.tol) == 1, eta.tol, sample(eta.tol, 1))
      eps <- CVaR.tol[eta.tol == eta] - min(CVaR.tol)
      lambda <- NA
    } else {
      eta <- HMCR.root(z, prob, k_star, ALPHA, tol.eta)
      lambda <- with(data, mean(
        surg(+f)/2 * pmax(+C - eta, 0)^k_star + surg(-f)/2 * pmax(-C - eta, 0)^k_star )^(1/k_star)
      )
      eps <- 0
    }
    const <- ifelse(k == Inf, CONST/2, CONST/(2 * k_star * lambda^{k_star - 1}))
    data$ZP <- with(data, const * pmax(+C - eta, 0)^k_star)
    data$ZN <- with(data, const * pmax(-C - eta, 0)^k_star)
    data$S <- with(data, ZP * surgN.gd(f) - ZN * surgN.gd(-f))

    if(verbose) {
      print(paste0("The ", iter, "-th iter: eta=", eta,
                   "; l2-norm of the slope=", slope.l2,
                   "; obj=", obj,
                   "; l2-norm of obj.gd=", obj.gd.l2))
    }

    if(obj.gd.l2 < tol.gd * max(slope.l2,1)) {
      message("Majorant update: gradient converges!"); break
    }
    # if(iter >= 2) {
    #   if(abs(obj - obj_old) < tol.obj * abs(obj_old)) {
    #     message("Majorant update: objective converges!"); break
    #   }
    # }

    par_old <- par
    obj_old <- obj
    L_old <- L(par_old, pen.l2, data)
    if(method.mm == "proj_gd") {
      for(iter.mm in 1:maxit.mm) {
        par0 <- par
        gd0 <- L.gd(par0, pen.l2, data)
        if(sum(gd0^2) < tol.gd.mm^2 * max(sum(par0^2),1)) break
        L0 <- L(par0, pen.l2, data)
        f0 <- with(data, par0[1] + as.vector(X %*% par0[-1]))
        W0 <- with(data, ZP * (f0 >= 0 & f0 <= 1) + ZN * (f0 >= -1 & f0 <= 0))
        Lips.lower <- 2/n * sv.X1^2 * mean(W0) + pen.l2
        Lips.upper <- 2/n * sv.X1^2 * max(W0) + pen.l2
        ### line-search
        Lips.ls <- Lips.lower
        repeat {
          par1 <- par0 - gd0/Lips.ls  # one-step gradient descent
          par <- c(par1[1], proj.l2(par1[-1], cnst.l2))
          if(L(par, pen.l2, data) <= L0 + sum(gd0 * (par - par0)) + Lips.ls/2 * sum((par - par0)^2)) {
            break  # Armijo condition: successful majorization
          } else if(Lips.ls >= Lips.upper) {
            warning(paste0("The ", iter, "-th majorant minimization: gradient descent failed to majorize!"))
            break
          }
          Lips.ls <- min(2 * Lips.ls, Lips.upper)
        }
        if(sum((par - par0)^2) <= tol.par.mm^2 * max(sum(par0^2),1)) break
      }
      if(iter.mm == maxit.mm) {
        warning(paste0("The ", iter, "-th majorant minimization: maximal iteration is attained!"))
      } else if(verbose) {
        print(paste0("Equivalent l2-penalty parameter=",
                     pen.l2_equiv <- max(sqrt(sum(L.gd(par,0,data)[-1]^2) / sum(par[-1]^2)), pen.l2)))
      }
    } else if(method.mm == "lbfgs") {
      pen.l2_equiv <- pen.l2
      for(ftol in 1/2 * 10^seq(0, -4, -1)) {
        optim <- lbfgs(
          L, L.gd, par, data = data, epsilon = tol.gd.mm, delta = tol.obj.mm,
          linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING", ftol = ftol, past = 20,
          max_iterations = maxit.mm, invisible = 1)
        if(optim$convergence != -998) break
      }
      if(optim$convergence == -998) {
        warning(paste0("The ", iter, "-th majorant minimization: L-BFGS maximum number of line-searches attained!"))
      } else if(ftol != 1/2) {
        warning(paste0("The ", iter, "-th majorant minimization: line-search accuracy parameter ", ftol))
      }
      if(optim$convergence == -997) {
        warning(paste0("The ", iter, "-th majorant minimization: L-BFGS doesn't converge!"))
      }
      par <- optim$par
    } else {
      stop('method.mm = "proj_gd" or "lbfgs"')
    }
    if(L_old - L(par, pen.l2, data) < eps) {
      par <- par_old
    } else if(sum((par - par_old)^2) < tol.par^2 * max(sum(par_old^2),1)) {
      message("Majorant update: parameter converges!"); break
    }
  }
  if(iter == maxit) {
    warning("Majorant update: maximal iteration attained!")
  }
  return(list(par = par, eta = eta, lambda = lambda,
              pen.l2 = ifelse(exists("pen.l2_equiv"), pen.l2_equiv, NA)))
}
DRITR <- function(
  data, par.init = NULL, cnst.l2 = Inf, pen.l2 = 1e-8, k = 2, CONST = 1, rho = NULL, method.Inf = "BSUM",
  method.mm = "proj_gd", tol.gd.mm = tol.gd, tol.obj.mm = tol.obj, tol.par.mm = tol.par, maxit.mm = maxit,
  tol.gd = 1e-5, tol.obj = 1e-5, tol.par = 1e-5, maxit = 200, tol.eta = 1e-5, tol.CVaR = 1e-8, verbose = F)
  # data = list of data (X, C)
  # par[1] = intercept
  # X doesn't include the all-ones (intercept) column
  # method.mm = "proj_gd" or "lbfgs"; cnst.l2 not available for lbfgs
{
  if(k != Inf & !is.null(rho)) CONST <- (k*(k-1)*rho + 1)^(1/k)
  n <- nrow(data$X)
  p <- ncol(data$X)
  if(is.null(par.init)) {
    par.init[2:(p+1)] <- rnorm(p) %>% { ifelse(cnst.l2 == Inf, 1, cnst.l2) * ./sqrt(sum(.^2)) }
    par.init[1] <- - sum(data$X[sample.int(n,1),] * par.init[-1])
  } else if(length(par.init) != p + 1) {
    stop(paste0("The length of an initial parameter should be ", p + 1))
  }
  par.init <- c(par.init[1], proj.l2(par.init[-1], cnst.l2))

  if(CONST < 1) {
    stop("CONST should be >= 1 or rho should be >= 0")
  } else if(CONST == 1) {
    message("Call DCA")
    return(DCA(data, par.init, cnst.l2, pen.l2,
               method.mm, tol.gd.mm, tol.obj.mm, tol.par.mm, maxit.mm,
               tol.gd, tol.obj, tol.par, maxit, verbose = F))
  } else if(k == Inf & method.Inf == "PrDCA") {
    message("Call PrDCA")
    return(PrDCA(data, par.init, cnst.l2, pen.l2, CONST,
                 method.mm, tol.gd.mm, tol.obj.mm, tol.par.mm, maxit.mm,
                 tol.gd, tol.obj, tol.par, maxit, tol.CVaR, verbose = F))
  } else if(k == Inf & method.Inf == "BSUM" | k > 1) {
    message("Call BSUM")
    return(BSUM(data, par.init, cnst.l2, pen.l2, k, CONST,
                method.mm, tol.gd.mm, tol.obj.mm, tol.par.mm, maxit.mm,
                tol.gd, tol.obj, tol.par, maxit, tol.eta, tol.CVaR, verbose = F))
  } else {
    if(!(method.Inf %in% c("BSUM", "PrDCA")))
      stop("method.Inf should be BSUM or PrDCA")
    else if(k <= 1)
      stop("k should be > 1 and <= Inf")
  }
}


