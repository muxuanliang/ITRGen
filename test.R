sim1 <- function(mu1=0, mu2=1.958){
  print(system.time(res <- foreach(index=1:500, .packages=c('ITRGen', 'grf', 'foreach'), .combine=rbind, .errorhandling='remove') %dopar% {
    set.seed(index)

    sample.size.train <- 1000
    sample.size.test <- 200
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

    mu.train <- 0 * rep(1, times=p)
    V.train <- V(p, 0)
    x.train <- mgcv::rmvn(sample.size.train, mu.train, V.train)
    tr.train <- rbinom(sample.size.train,1,0.5)
    y.train <- (tr.train-0.5) * (x.train[,2]-((x.train[,1])^3-2*x.train[,1]))+1+apply(x.train,1,mean)+rnorm(sample.size.train)

    mu.test <- c(mu1, mu2, rep(0, times=p-2))
    V.test <- V(p, 0)
    x.test <- mgcv::rmvn(sample.size.test, mu.test, V.test)

    # estimate density ratio
    #fit_weighted <- ITRFitAll(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), propensity = rep(0.5, times=sample.size.train), is.weight = TRUE, x.test = x.test, test = FALSE)
    #fit_null <- ITRFitAll(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), propensity = rep(0.5, times=sample.size.train), test = FALSE)
    fit_w <- ContrastITR(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), is.weight = TRUE, x.test = x.test)
    fit_c <- ContrastITR(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), is.weight = FALSE, x.test = x.test)
    # one directional
    #beta <- (fit_null$fit[[1]]$fit$beta[,fit_null$fit[[1]]$fit$lambda==fit_null$fit[[1]]$fit$lambda.min]+fit_null$fit[[2]]$fit$beta[,fit_null$fit[[2]]$fit$lambda==fit_null$fit[[2]]$fit$lambda.min])/2
    #density.ratio.one <- densratio::densratio(x.train %*% beta, x.test %*% beta)
    #w.weight.one <- pmin(1/density.ratio.one$compute_density_ratio(x.train %*% beta), rep(bdd, times=sample.size.train))
    #fit_weighted.one <- ITRFitAll(data=list(predictor = x.train, treatment = tr.train, outcome=y.train), propensity = rep(0.5, times=sample.size.train), sample.weight = w.weight.one)

    # test dataset
    x.test <- mgcv::rmvn(10^6, mu.test, V.test)
    #d_weighted.one <- sign(predict(fit_weighted.one$fit[[1]]$fit, newx = x.test, s=fit_weighted.one$fit[[1]]$fit$lambda.min)+predict(fit_weighted.one$fit[[2]]$fit, newx = x.test, s=fit_weighted.one$fit[[2]]$fit$lambda.min))
    #d_weighted <- sign(predict(fit_weighted$fit[[1]]$fit, newx = x.test, s=fit_weighted$fit[[1]]$fit$lambda.min)+predict(fit_weighted$fit[[2]]$fit, newx = x.test, s=fit_weighted$fit[[2]]$fit$lambda.min))
    #d_null <- sign(predict(fit_null$fit[[1]]$fit, newx = x.test, s=fit_null$fit[[1]]$fit$lambda.min)+predict(fit_null$fit[[2]]$fit, newx = x.test, s=fit_null$fit[[2]]$fit$lambda.min))
    d_s <- sign(fit_c$standard_fit[1]+x.test%*%fit_c$standard_fit[-1])
    d_w <- sign(fit_w$standard_fit[1]+x.test%*%fit_w$standard_fit[-1])
    d_c <- sign(fit_c$dr_fit[1]+x.test%*%fit_c$dr_fit[-1])
    #value_weighted.one <- mean((x.test[,2]-((x.test[,1])^3-2*x.test[,1])) * d_weighted.one)
    value_s <- mean((x.test[,2]-((x.test[,1])^3-2*x.test[,1])) * d_s)
    value_w <- mean((x.test[,2]-((x.test[,1])^3-2*x.test[,1])) * d_w)
    value_c <- mean((x.test[,2]-((x.test[,1])^3-2*x.test[,1])) * d_c)

    c(value_s, value_w, value_c)
  }))
  save(res, file = paste0("/mnt/c/Users/lmx19/Documents/Simulations/ITRGen/case1_", mu1, "_", mu2,"_", sample.size.test,".RData"))
  apply(res,2,mean)
  apply(res,2,sd)
}

library(doParallel)
n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

mu1_seq <- mu2_seq <- c(1.958, 1.469, 0.734, 0)#seq(-2.448, 2.448, length.out = 21)

for (mu1 in mu1_seq){
  for (mu2 in mu2_seq){
    if (mu1==0){
      sim1(mu1 = mu1, mu2 = mu2)
    }
    if ((mu1==1.958) & (mu2 %in% c(1.469, 1.958))){
      sim1(mu1 = mu1, mu2 = mu2)
    }
  }
}


stopCluster(cl)





