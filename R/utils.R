# loss gets the loss of the chosen loss
loss <- function(x, loss_type){
  switch(loss_type,
         logistic = log(1+exp(-x)),
         exponential = exp(-x),
         zeroOne = (x<0)
  )
}

# derivative gets the derivative of the chosen loss
derivative <- function(x, loss_type){
  switch(loss_type,
         logistic = -exp(-x)/(1+exp(-x)),
         exponential = -exp(-x)
  )
}

# hessian gets the hessian of the chosen loss
hessian <- function(x, loss_type){
  switch(loss_type,
         logistic = exp(-x)/(1+exp(-x))^2,
         exponential = exp(-x)
  )
}

# ks gets the kernel estimation
ks <- function(xx, yy, xx.test){
  nobs <- nrow(as.matrix(xx))
  nvars <- ncol(as.matrix(xx))
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
  wm <- function(t){
    if (ncol(as.matrix(xx))==1){
      weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
    }
    weighted.mean(yy, weight)
  }
  if (nrow(as.matrix(xx.test))==1) {
    yy.test <- wm(xx.test)
  } else {
    if (ncol((as.matrix(xx.test)))==1){
      yy.test <- sapply(as.matrix(xx.test), function(t){
        wm(t)
      })
    } else {
      yy.test <- apply(as.matrix(xx.test),1,function(t){
        wm(t)
      })
    }
  }
  yy.test
}

# model.gam gets the gam regression function given data
model.gam <- function(data){
  p <- dim(data$predictor)[2]
  expr <- "mgcv::gam(outcome~"
  for (i in 1:(p-1)){
    expr <- paste0(expr, "s(predictor[,",i,"])+")
  }
  expr <- paste0(expr, "s(predictor[,",p,"]), data = data,method ='REML')")
  expr
}

# getOutcomeModel contains the outcome regression model
getOutcomeModel <- function(data, method=c('lm', 'glmnet', 'kernel', 'others'), sampleSplitIndex, Formula = NULL, predictAll = FALSE, screeningMethod="SIRS", outcomeScreeningFamily='Gaussian'){
  p <- dim(data$predictor)[2]
  size <- dim(data$predictor)[1]
  fit <- NULL
  supp <- NULL
  dataPredict <- NULL
  dataPredict$control=data$predictor[sampleSplitIndex,]
  dataPredict$treatment=data$predictor[sampleSplitIndex,]
  if (predictAll){
    dataPredict$control=data$predictor
    dataPredict$treatment=data$predictor
  }
  prediction <- NULL
  dataControl <- list(predictor=data$predictor[(!sampleSplitIndex) & (data$treatment==FALSE),], outcome=data$outcome[(!sampleSplitIndex) & (data$treatment==FALSE)])
  dataTreatment <- list(predictor=data$predictor[(!sampleSplitIndex) & (data$treatment==TRUE),], outcome=data$outcome[(!sampleSplitIndex) & (data$treatment==TRUE)])
  if (0.05*size >= p){
    supp$control <- supp$treatment <- rep(TRUE, times = p)
    if ((method != 'lm')&&(method != 'glmnet')){
      ans1 <- screening(data$predictor[data$treatment==FALSE,], data$outcome[data$treatment==FALSE], method = screeningMethod, family = outcomeScreeningFamily)
      ans2 <- screening(data$predictor[data$treatment==TRUE,], data$outcome[data$treatment==TRUE], method = screeningMethod, family = outcomeScreeningFamily)
      if (screeningMethod == 'glmnet'){
        supp$control <- ans1
        supp$treatment <- ans2
      } else {
        supp$control <- supp$treatment <- (ans1 <= floor(p/2))|(ans2 <= floor(p/2))
      }
      dataControl$predictor <- dataControl$predictor[,supp$control]
      dataTreatment$predictor <- dataTreatment$predictor[,supp$treatment]
    }
  }
  if ((0.05*size < p) || (method == 'glmnet')) {
    fit$control <- glmnet::cv.glmnet(x = dataControl$predictor, y = dataControl$outcome)
    fit$treatment <- glmnet::cv.glmnet(x = dataTreatment$predictor, y = dataTreatment$outcome)
    supp$control <- abs(fit$control$glmnet.fit$beta[,fit$control$glmnet.fit$lambda==fit$control$lambda.min])>0
    supp$treatment <- abs(fit$treatment$glmnet.fit$beta[,fit$treatment$glmnet.fit$lambda==fit$treatment$lambda.min])>0
    if ((method != 'lm')&&(method != 'glmnet')){
      ans1 <- screening(data$predictor[data$treatment==FALSE,], data$outcome[data$treatment==FALSE], method = screeningMethod, family = outcomeScreeningFamily)
      ans2 <- screening(data$predictor[data$treatment==TRUE,], data$outcome[data$treatment==TRUE], method = screeningMethod, family = outcomeScreeningFamily)
      if (screeningMethod == 'glmnet'){
        supp$control <- ans1
        supp$treatment <- ans2
      } else {
        supp$control <- supp$treatment <- (ans1 <= 5)|(ans2 <= 5)
      }
    }
    dataControl$predictor <- dataControl$predictor[,supp$control]
    dataTreatment$predictor <- dataTreatment$predictor[,supp$treatment]
    prediction$control <- predict(fit$control, newx = dataPredict$control, s=fit$control$lambda.min)
    prediction$treatment <- predict(fit$treatment, newx = dataPredict$treatment, s=fit$treatment$lambda.min)
  }

  if (is.null(Formula)){
    Formula <- function(support){
      expr <- (outcome ~ predictor)
      if (sum(support)==0){
        expr <- (outcome ~ 1)
      }
    expr
    }
  }
  fit <- NULL
  dataPredict <- NULL
  dataPredict$control=data$predictor[sampleSplitIndex,supp$control]
  dataPredict$treatment=data$predictor[sampleSplitIndex,supp$treatment]
  if (predictAll){
    dataPredict$control=data$predictor[,supp$control]
    dataPredict$treatment=data$predictor[,supp$treatment]
  }
  if ((method == 'lm')||(method == 'glmnet')){
    if (sum(supp$control) > 0){
      fit$control <- lm(Formula(supp$control), data = dataControl)
      prediction$control <- predict(fit$control, newdata = list(predictor=dataPredict$control))
    }
    if (sum(supp$treatment) > 0){
      fit$treatment <- lm(Formula(supp$treatment), data = dataTreatment)
      prediction$treatment <- predict(fit$treatment, newdata = list(predictor=dataPredict$treatment))
    }
  } else if (method == 'kernel') {
    if (sum(supp$control) > 0){
    prediction$control <- ks(dataControl$predictor, dataControl$outcome, dataPredict$control)
    }
    if (sum(supp$treatment) > 0){
    prediction$treatment <- ks(dataTreatment$predictor, dataTreatment$outcome, dataPredict$treatment)
    }
  } else {
    if (sum(supp$control) > 0){
    fit$control <- eval(parse(text=model.gam(dataControl)))
    prediction$control <- predict(fit$control, newdata = list(predictor=dataPredict$control))
    }
    if (sum(supp$treatment) > 0){
    fit$treatment <- eval(parse(text=model.gam(dataTreatment)))
    prediction$treatment <- predict(fit$treatment, newdata = list(predictor=dataPredict$treatment))
    }
  }
  prediction
}

# getPropensityModel contains the outcome regression model
getPropensityModel <- function(data, method=c('lm', 'glmnet', 'kernel'), sampleSplitIndex, Formula = NULL, predictAll = FALSE, screeningMethod="SIRS"){
  p <- dim(data$predictor)[2]
  size <- dim(data$predictor)[1]
  fit <- NULL
  supp <- NULL
  dataPredict <- NULL
  dataPredict=data$predictor[sampleSplitIndex,]
  if (predictAll){
    dataPredict=data$predictor
  }
  prediction <- NULL
  dataTrain <- list(predictor=data$predictor[(!sampleSplitIndex),], treatment=data$treatment[(!sampleSplitIndex)])
  if (0.05*size >= p){
    supp$control <- supp$treatment <- rep(TRUE, times = p)
    if ((method != 'lm')&&(method != 'glmnet')){
      ans <- screening(data$predictor, data$treatment, method = screeningMethod, family = 'binomial')
      if (screeningMethod == 'glmnet'){
        supp <- ans
      } else {
        supp <- (ans <= p/2)
      }
    }
  }
  if ((0.05*size < p) || (method == 'glmnet')) {
    fit <- glmnet::cv.glmnet(x = dataTrain$predictor, y = dataTrain$treatment, family='binomial')
    supp <- abs(fit$glmnet.fit$beta[,fit$glmnet.fit$lambda==fit$lambda.min])>0
    if ((method != 'lm')&&(method != 'glmnet')){
      ans <- screening(data$predictor, data$treatment, method = screeningMethod, family = 'binomial')
      if (screeningMethod == 'glmnet'){
        supp <- ans
      } else {
        supp <- (ans <= 5)
      }
    }
    dataTrain$predictor <- dataTrain$predictor[,supp]
    prediction <- predict(fit, newx = dataPredict, type='response', s=fit$lambda.min)
  }

  if (is.null(Formula)){
    Formula <- function(support){
      expr <- (treatment ~ predictor)
      if (sum(support)==0){
        expr <- (treatment ~ 1)
      }
      expr
    }
  }
  fit <- NULL
  dataPredict <- NULL
  dataPredict=data$predictor[sampleSplitIndex,supp]
  if (predictAll){
    dataPredict=data$predictor[,supp]
  }
  if ((method == 'lm')||(method == 'glmnet')){
    if (sum(supp) > 0){
      fit <- glm(Formula(supp), family=binomial, data = dataTrain)
      prediction <- predict(fit, newdata = list(predictor=dataPredict), type="response")
    }
  } else if (method == 'kernel') {
    if (sum(supp) > 0){
      prediction <- ks(dataTrain$predictor, dataTrain$treatment, dataPredict)
      prediction <- (prediction > 0.9) * 0.9 + (prediction < 0.1) * 0.1 + (prediction < 0.9) * (prediction > 0.1) * prediction
    }
  }
  prediction
}

# screening
screening <- function(x, y, method='glmnet', family='Gaussian'){
  var <- apply(x, 2, sd)
  supp <- order(var, decreasing = TRUE)
  if (method=='glmnet'){
    fit <- glmnet::cv.glmnet(x, y, family = family)
    coef <- fit$glmnet.fit$beta[,fit$lambda==fit$lambda.min]
    supp <- (abs(coef)>0)
  } else {
    fit <- VariableScreening::screenIID(x, y, method=method)
    supp <- fit$rank
  }
  supp
}

library(dplyr)

##### True objectives
OBJ <- function(par, data, k = 2, CONST = 1, rho = NULL)
  # objective
  # data = list of data (X, C)
  # par[1] = intercept
  # X doesn't include the all-ones (intercept) column
{
  f <- with(data, par[1] + as.vector(X %*% par[-1]))
  if(k != Inf & !is.null(rho)) CONST <- (k*(k-1)*rho + 1)^(1/k)

  if(CONST < 1) {
    stop("CONST should be >= 1 and rho should be >= 0")
  } else if(CONST == 1) {  # ERM
    return(with(data, mean(C * (surg(f) - 1))))
  }
  ALPHA = 1 - 1/CONST
  k_star <- k/(k-1)
  z <- with(data, c(C, -C))
  prob <- c(surg(f), surg(-f))

  if(k == Inf) {
    obj <- min(CVaR(,z, prob, ALPHA))
  } else {
    eta <- HMCR.root(z, prob, k_star, ALPHA)
    obj <- HMCR(eta, z, prob, k_star, ALPHA)
  }
  return(obj)
}
OBJ.gd <- function(par, data, k = 2, CONST = 1, rho = NULL)
  # objective gradient: colMeans(((ZP - ZN) * surg.gd(f)) * cbind(1,X))
  # data = list of data (X, C)
  # par[1] = intercept
  # X doesn't include the all-ones (intercept) column
{
  f <- with(data, par[1] + as.vector(X %*% par[-1]))
  if(k != Inf & !is.null(rho)) CONST <- (k*(k-1)*rho + 1)^(1/k)

  if(CONST < 1) {
    stop("CONST should be >= 1 and rho should be >= 0")
  } else if(CONST == 1) {  # ERM
    return(with(data, colMeans(C * surg.gd(f) * cbind(1,X))))
  }
  ALPHA <- 1 - 1/CONST
  k_star <- ifelse(k == Inf, 1, k/(k-1))
  z <- with(data, c(C, -C))
  prob <- c(surg(f), surg(-f))

  if(k == Inf) {
    eta <- VaR(z, prob, ALPHA)[1]
    lambda <- 1
  } else {
    eta <- HMCR.root(z, prob, k_star, ALPHA)
    lambda <- with(data, mean( surg(+f)/2 * pmax(+C - eta, 0)^k_star +
                                 surg(-f)/2 * pmax(-C - eta, 0)^k_star )^(1/k_star))
  }
  const <- CONST / (2 * k_star * lambda^{k_star - 1})
  data$ZP <- with(data, const * pmax(+C - eta, 0)^k_star)
  data$ZN <- with(data, const * pmax(-C - eta, 0)^k_star)
  return(with(data, colMeans((ZP - ZN) * surg.gd(f) * cbind(1,X))))
}
##### ERM majorant objectives
L_ERM <- function(par, pen.l2 = 0, data)
  # majorant objective: mean(CP * surgP(f) + CN * surgN(f) - S * f) + pen
  # data = list of data (X, C, S)
  # par[1] = intercept
  # X doesn't include the all-ones (intercept) column
{
  f   <- with(data, par[1] + as.vector(X %*% par[-1]))
  pen <- pen.l2/2 * sum(par[-1]^2)
  return(with(data, mean(pmax(C, 0) * surgP(f) + pmax(-C, 0) * surgN(f) - S * f)) + pen)
}
L_ERM.gd <- function(par, pen.l2 = 0, data)
  # majorant gradient: colMeans((pmax(C, 0) * surgP.gd(f) + pmax(-C, 0) * surgN.gd(f) - S) * cbind(1,X))) + pen.gd
  # data = list of data (X, C, S)
  # par[1] = intercept
  # X doesn't include the all-ones (intercept) column
{
  f <- with(data, par[1] + as.vector(X %*% par[-1]))
  pen.gd <- pen.l2 * c(0, par[-1])
  return(with(data, colMeans((pmax(C, 0) * surgP.gd(f) + pmax(-C, 0) * surgN.gd(f) - S) * cbind(1,X))) + pen.gd)
}
##### BSUM majorant objectives
L <- function(par, pen.l2 = 0, data)
  # majorant objective: mean(ZP * surgP(f) + ZN * surgP(-f) - S * f) + pen
  # data = list of data (X, ZP, ZN, S)
  # par[1] = intercept
  # X doesn't include the all-ones (intercept) column
{
  f   <- with(data, par[1] + as.vector(X %*% par[-1]))
  pen <- pen.l2/2 * sum(par[-1]^2)
  return(with(data, mean(ZP * surgP(f) + ZN * surgP(-f) - S * f)) + pen)
}
L.gd <- function(par, pen.l2 = 0, data)
  # majorant gradient: colMeans((ZP * surgP.gd(f) - ZN * surgP.gd(-f) - S) * cbind(1,X))) + pen.gd
  # data = list of data (X, ZP, ZN, S)
  # par[1] = intercept
  # X doesn't include the all-ones (intercept) column
{
  f       <- with(data, par[1] + as.vector(X %*% par[-1]))
  pen.gd  <- pen.l2 * c(0, par[-1])
  return(with(data, colMeans((ZP * surgP.gd(f) - ZN * surgP.gd(-f) - S) * cbind(1,X))) + pen.gd)
}

library(dplyr)
##### Smoothed ramp losses
surgP <- function(u) { ifelse(u <= 1, ifelse(u <=  0,  1-2*u, (1-u)^2), 0) }
surgN <- function(u) { ifelse(u <= 0, ifelse(u <= -1, -1-2*u,     u^2), 0) }
surg  <- function(u) { surgP(u) - surgN(u) }
surgP.gd <- function(u) { ifelse(u <= 1, ifelse(u <=  0, -2, -2+2*u), 0) }
surgN.gd <- function(u) { ifelse(u <= 0, ifelse(u <= -1, -2,    2*u), 0) }
surg.gd  <- function(u) { surgP.gd(u) - surgN.gd(u) }
##### Coherent risk measures
HMCR <- function(eta=NULL, z, prob=rep(1/length(z),length(z)), order=2, alpha=0) {
  prob <- prob/sum(prob)
  if(is.null(eta)) { eta <- z }
  sapply(eta, function(eta) eta + 1/(1-alpha) * sum(prob * pmax(z - eta, 0)^order)^(1/order))
}
HMCR.gd <- function(eta=NULL, z, prob=rep(1/length(z),length(z)), order=2, alpha=0) {
  prob <- prob/sum(prob)
  if(is.null(eta)) { eta <- z }
  sapply(eta, function(eta)
    if(all(z <= eta)) return(Inf) else return(
      1 - 1/(1-alpha) * sum(prob * pmax(z - eta, 0)^order)^(1/order-1) *
        sum(prob * pmax(z - eta, 0)^(order-1))
    ))
}
HMCR.root <- function(z, prob=rep(1/length(z),length(z)), order=2, alpha=0, tol.eta=1e-5) {
  prob <- prob/sum(prob)
  sapply(alpha, function(alpha)
    if(alpha == 0) return(-Inf) else return(
      uniroot(HMCR.gd, z = z, prob = prob, order = order, alpha = alpha,
              interval = c((min(z)-(1-alpha)*max(z))/alpha, max(z*(prob>0))-tol.eta), tol = tol.eta)$root
    ))
}
CVaR <- function(eta=NULL, z, prob=rep(1/length(z),length(z)), alpha=0) {
  prob <- prob/sum(prob)
  if(is.null(eta)) {
    if(length(z) == 1) return(z) else {
      id.sort <- sort(z, decreasing = T, index.return = T)$ix
      CVaR_sort <- z[id.sort[1]]
      CCDF <- 0  # complementary CDF
      for(i in 2:length(z)) {
        z_sort_diff <- z[id.sort[i-1]] - z[id.sort[i]]
        CCDF <- CCDF + prob[id.sort[i-1]]
        CVaR_sort[i] <- CVaR_sort[i-1] - (1 - CCDF/(1-alpha)) * z_sort_diff
      }
      CVaR <- numeric()
      CVaR[id.sort] <- CVaR_sort
      return(CVaR)
    }
  } else return(
    sapply(eta, function(eta) eta + 1/(1-alpha) * sum(prob * pmax(z - eta,0)))
  )
}
VaR <- function(z, prob=rep(1/length(z),length(z)), alpha=0, tol.CVaR=0) {
  prob <- prob/sum(prob)
  CVaR(,z, prob, alpha) %>% { z[which(. <= min(.) + tol.CVaR)] }
}
##### Projection on the l2-ball
proj.l2 <- function(x, cnst.l2 = Inf) {
  return(x / pmax(sqrt(sum(x^2))/cnst.l2, 1))
}

