sim2 <- function(xi=0.1){
  print(system.time(res <- foreach(index=1:500, .packages=c('ITRGen', 'grf', 'foreach'), .combine=rbind, .errorhandling='remove') %dopar% {
    set.seed(index)

    sample.size.train <- 1000
    sample.size.test <- 50
    p <- 10

    V <- function(p, rate = 0.5){
      V.matrix <- array(0, c(p,p))
      for (i in 1:p){
        for (j in 1:p){
          V.matrix[i,j] <- rate ^ (abs(i-j))
        }
      }
      V.matrix
    }

    mu1 <- 0.5 * c(-1, 1, rep(0, times=p-2))
    mu0 <- (- mu1)
    xi.train <- rbinom(sample.size.train, 1, 0.75)
    V.train <- V(p, 0)
    x.train <- matrix(rnorm(p * sample.size.train), ncol = p) + xi.train %o% mu1 + (1-xi.train) %o% mu0
    tr.train <- rbinom(sample.size.train,1,0.5)
    y.train <- (tr.train-0.5) * ((-1.5) * (2*xi.train -1)-2 * x.train[,1]+x.train[,2])+1+apply(x.train,1,mean)+rnorm(sample.size.train)

    V.test <- V(p, 0)
    xi.test <- rbinom(sample.size.test, 1, xi)
    x.test <- matrix(rnorm(p * sample.size.test), ncol = p) + xi.test %o% mu1 + (1-xi.test) %o% mu0

    # estimate density ratio
    #density.ratio <- densratio::densratio(x.train, x.test)
    #w.weight <- pmin(1/density.ratio$compute_density_ratio(x.train), rep(bdd, times=sample.size.train))
    #fit_weighted <- ITRFitAll(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), propensity = rep(0.5, times=sample.size.train), is.weight = TRUE, x.test = x.test)
    #fit_null <- ITRFitAll(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), propensity = rep(0.5, times=sample.size.train))
    fit_w <- ContrastITR(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), is.weight = TRUE, x.test = x.test)
    fit_c <- ContrastITR(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), is.weight = FALSE, x.test = x.test)
    # one directional
    #beta <- (fit_null$fit[[1]]$fit$beta[,fit_null$fit[[1]]$fit$lambda==fit_null$fit[[1]]$fit$lambda.min]+fit_null$fit[[2]]$fit$beta[,fit_null$fit[[2]]$fit$lambda==fit_null$fit[[2]]$fit$lambda.min])/2
    #density.ratio.one <- densratio::densratio(x.train %*% beta, x.test %*% beta)
    #w.weight.one <- pmin(1/density.ratio.one$compute_density_ratio(x.train %*% beta), rep(bdd, times=sample.size.train))
    #fit_weighted.one <- ITRFitAll(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), propensity = rep(0.5, times=sample.size.train), sample.weight = w.weight.one)

    # test dataset
    xi.test <- rbinom(10^6, 1, xi)
    x.test <- matrix(rnorm(p * 10^6), ncol = p) + xi.test %o% mu1 + (1-xi.test) %o% mu0
    #d_weighted.one <- sign(predict(fit_weighted.one$fit[[1]]$fit, newx = x.test, s=fit_weighted.one$fit[[1]]$fit$lambda.min)+predict(fit_weighted.one$fit[[2]]$fit, newx = x.test, s=fit_weighted.one$fit[[2]]$fit$lambda.min))
    #d_weighted <- sign(predict(fit_weighted$fit[[1]]$fit, newx = x.test, s=fit_weighted$fit[[1]]$fit$lambda.min)+predict(fit_weighted$fit[[2]]$fit, newx = x.test, s=fit_weighted$fit[[2]]$fit$lambda.min))
    #d_null <- sign(predict(fit_null$fit[[1]]$fit, newx = x.test, s=fit_null$fit[[1]]$fit$lambda.min)+predict(fit_null$fit[[2]]$fit, newx = x.test, s=fit_null$fit[[2]]$fit$lambda.min))
    d_s <- sign(fit_c$standard_fit[1]+x.test%*%fit_c$standard_fit[-1])
    d_w <- sign(fit_w$standard_fit[1]+x.test%*%fit_w$standard_fit[-1])
    d_c <- sign(fit_c$dr_fit[1]+x.test%*%fit_c$dr_fit[-1])
    #value_weighted.one <- mean(((-1.5) * (2*xi.test -1)-2 * x.test[,1]+x.test[,2]) * d_weighted.one)
    value_weighted <- mean(((-1.5) * (2*xi.test -1)-2 * x.test[,1]+x.test[,2]) * d_weighted)
    value_null <- mean(((-1.5) * (2*xi.test -1)-2 * x.test[,1]+x.test[,2]) * d_null)
    value_c <- mean(((-1.5) * (2*xi.test -1)-2 * x.test[,1]+x.test[,2]) * d_c)

    c(value_weighted, value_null, value_c)
  }))
  save(res, file = paste0("/mnt/c/Users/lmx19/Documents/Simulations/ITRGen/case2_", xi,".RData"))
  #save(res, file = paste0("~/Simulations/ITRGen/case2_", xi,".RData"))
  apply(res,2,mean)
  apply(res,2,sd)
}

library(doParallel)
n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

xi_seq <- c(0.1, 0.25, 0.5, 0.75, 0.9)#seq(-2.448, 2.448, length.out = 21)

for (xi in xi_seq){
    sim2(xi=xi)
}

stopCluster(cl)





