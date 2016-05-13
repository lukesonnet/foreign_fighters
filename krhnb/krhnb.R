#----------
# Authors: George Derpanopoulos and Luke Sonnet
#
# Project: Kernel Regularized Hurdle Negative Binomial
#
# Purpose: Provide all of the functions needed to fit and analyze KRHNB. Also
#          includes an example using the imputed Foreign Fighters data and more
#----------

library(pdist)
library(pscl)
#library(foreach)
#todo: add parallelization for bootstrap
library(ggplot2)
library(reshape2)
#todo: replace the cross validation summing each time to just a vector of yhats and then compute overall prediction error

#----------
# Functions
#----------

## This is the main function. Fits the model by preparing the data, searching
## for lambdas through CV, minimizing the target function, and returning an
## object that contains all of the information about the fitting that is
## necessary for future analysis. I borrow liberally from the KRLS package
krhnb <- function(X = NULL,
                  y = NULL,
                  whichkernel = "gaussian",  
                  folds = 5,
                  lambda1 = NULL,
                  lambda2 = NULL,
                  lambdafolds = 5,
                  lambdastart = c(0, 0),
                  lambda1s = NULL,
                  lambda2s = NULL,
                  sigma = NULL,
                  con = list(maxit=500),
                  printout = TRUE,
                  returnoptim = TRUE,
                  hessian = FALSE,
                  effects = TRUE,
                  bootstrap = TRUE,
                  bootstrapits = 100
                  #parallel = FALSE
                  ) {
  
  ## Prepare the data
  X <- as.matrix(X)
  y <- as.matrix(y)
  
  ## Input validation
  if (is.numeric(X)==FALSE) {
    stop("X must be numeric")
  }      
  if (is.numeric(y)==FALSE) {
    stop("y must be numeric")
  }
  if (sd(y)==0) {
    stop("y is a constant")        
  }  
  if (sum(is.na(X))>0){
    stop("X contains missing data")
  }
  if (sum(is.na(y))>0){
    stop("y contains missing data")
  }
  if (!all(y%%1==0)) {
    stop("y is not count data")
  }
  
  ## Default sigma to the number of features
  if(is.null(sigma)){
    sigma <- 2*ncol(X)
  } else{
    stopifnot(is.vector(sigma),
              length(sigma)==1,
              is.numeric(sigma),
              sigma>0)        
  }
  
  ## Scale data
  X.init <- X
  X.init.sd <- apply(X.init, 2, sd)
  X.init.mean <- apply(X.init, 2, mean)
  if (sum(X.init.sd == 0)){
    stop("at least one column in X is a constant, please remove the constant(s)")
  }
  X <- scale(X, center = X.init.mean, scale = X.init.sd)  
  
  ## Compute kernel matrix
  K <- NULL
  if(whichkernel=="gaussian"){ K <- gaussKernel(X, sigma = sigma)}
  if(whichkernel=="linear"){ K <- tcrossprod(X) }
  if(whichkernel=="poly2"){ K <- (tcrossprod(X)+1)^2 }
  if(whichkernel=="poly3"){ K <- (tcrossprod(X)+1)^3 }
  if(whichkernel=="poly4"){ K <- (tcrossprod(X)+1)^4 }
  if(is.null(K)){ stop("No valid Kernel specified") }
  
  ## Get starting value for theta
  hnb <- NULL
  try(hnb <- hurdle(V1 ~ ., data = as.data.frame(cbind(y, X)),
                dist = "negbin", y = FALSE)$theta)
  initTheta <- hnb
  
  ## Grid search for lambdas, this implies that if either value is not specified
  ## it will attempt a grid search. Better to use lambda1s and lambda2s even
  ## if you only want to look at one value of one of the vectors
  if (is.null(lambda1) | is.null(lambda2)) {
    chunks <- chunk(sample(nrow(K)), lambdafolds)
    
    if (is.null(lambda1s)) {
      lambda1s <- c(0.01, 0.05, 0.1, 0.5, 1, 10, 20, 50)
    } 
    if (is.null(lambda2s)) {
      lambda2s <- c(0.01, 0.05, 0.1, 0.5, 1, 10, 20, 50)
    }
    
    lambdaGrid <- cbind(expand.grid(lambda1s, lambda2s), NA)
    
    for(i in 1:nrow(lambdaGrid)){
      lambdaGrid[i, 3] <- lambda.fn(par = log(as.numeric(lambdaGrid[i,1:2])),
                                    X=X.init, y=y,
                                    folds=lambdafolds, chunks=chunks)
      print(lambdaGrid[i, ]) #todo: change printout options
    }
    print(which.min(lambdaGrid[, 3]))
    lambdas.min <- lambdaGrid[which.min(lambdaGrid[, 3]), 1:2]
    lambda1 <- as.numeric(lambdas.min[1])
    lambda2 <- as.numeric(lambdas.min[2])
    
  } else {
    # check user specified lambda
    stopifnot(is.vector(lambda1),
              length(lambda1)==1,
              is.numeric(lambda1),
              lambda1>0,
              is.vector(lambda2),
              length(lambda2)==1,
              is.numeric(lambda2),
              lambda2>0)  
  }

  ## Solve!
  ## Initialize parameters
  pars <- rep(0, 2*nrow(K) + 3)
  if (!is.null(initTheta)) pars[length(pars) - 2] <- log(initTheta)

  parfitted <- solveForC(pars = pars,
                         y=y, K=K, lambda1=lambda1, lambda2=lambda2, con=con,
                         hessian=hessian)
  
  chat1 <- parfitted$chat1
  chat2 <- parfitted$chat2
  thetahat <- parfitted$thetahat
  beta0_1 <- parfitted$beta0_1
  beta0_2 <- parfitted$beta0_2
  
  #print(initTheta)
  #print(thetahat)
  
  ## Fitted Ys
  yfitted <- predictY(c1=chat1, c2=chat2, theta=thetahat, beta0_1=beta0_1,
                      beta0_2=beta0_2, K=K)

  
  # Getting marginal effects
  if(effects){
    eff <- effects(X, parfitted)
    effectmat <- eff$effectDF
    onesmat <- eff$onesDF
    countmat <- eff$countDF
    meaneffect <- apply(effectmat, 2, mean)
  }

  ## Bootstrapping
  if(bootstrap){
    its <- bootstrapits
    bootmfx <- matrix(NA, nrow=its, ncol=ncol(X.init), dimnames=list(NULL, colnames(X.init)))
    bootonesmfx <- matrix(NA, nrow=its, ncol=ncol(X.init), dimnames=list(NULL, colnames(X.init)))
    bootcountmfx <- matrix(NA, nrow=its, ncol=ncol(X.init), dimnames=list(NULL, colnames(X.init)))
    bootmfx.median <- matrix(NA, nrow=its, ncol=ncol(X.init), dimnames=list(NULL, colnames(X.init)))
  
    for(i in 1:its){
      bootc <- NULL
      bootc$theta <- 1e-8
      ## The overdispersion parameter sometimes gets very small, resulting in the
      ## algorithm being stuck in a very bad local minima where predictions are
      ## far too high. Thus we resample until we avoid that local minima
      ## todo: allow users to set their own cutoff
      while(bootc$theta < 1e-2){
        ## get bootstrap sample
        samp <- sample(length(y), replace=T)
        
        ## use lambdas from full model, will limit variability
        newK <- gaussKernel(X = scale(X.init[samp, ]), sigma = ncol(X))
        
        bootc <- solveForC(y = y[samp],
                           K = newK,
                           lambda1=lambda1,
                           lambda2=lambda2)
            
        if(bootc$theta >= 1e-2){
          eff <- effects(X.init[samp, ], bootc)
          bootmfx[i, ] <- apply(eff$effectDF, 2, mean)
          bootmfx.median[i, ] <- apply(eff$effectDF, 2, median)
          bootonesmfx[i, ] <- apply(eff$onesDF, 2, mean)
          bootcountmfx[i, ] <- apply(eff$countDF, 2, mean)
          if (i %% 10 == 0) print(i)
        } else {
          print(c(i, bootc$theta))
        }
      }
    }
  }

  # return
  z <- list(K=K,
            chat1=chat1,
            chat2=chat2,
            thetahat=thetahat,
            beta0_1 = beta0_1,
            beta0_2 = beta0_2,
            fitted=yfitted,
            X=X.init,
            y=y,
            sigma=sigma,
            lambda1=lambda1,
            lambda2=lambda2,
            optim.obj = parfitted,
            kernel=whichkernel,
            effectmat = effectmat,
            onesmat = onesmat,
            countmat = countmat,
            meaneffect = meaneffect,
            bootmfx = if (bootstrap) bootmfx else NULL,
            bootmfx.median = if (bootstrap) bootmfx.median else NULL,
            bootonesmfx = if (bootstrap) bootonesmfx else NULL,
            bootcountmfx = if (bootstrap) bootcountmfx else NULL
            )
  
  class(z) <- "krhnb"  
  
  ## Output the bootstrapped median summary
  print(summary(z, type = "mean"))
  return(z)
}


## The Kernel Regularized Hurdle Negative Binomial target function to be minimized
## Parameters:
##   'par' - the parameters to be optimized, contains C and log(Theta)
##   'K' - the Kernel matrix
##   'y' - the outcome variable
##   'lambda' - the regularizing parameter
## Values:
##   'r' - The penalized log-likelihood given 'y' and 'K'
krhnb.fn <- function(par, K, y, lambda1=0.5, lambda2=0.5) {

  coef1  <- par[1:nrow(K)]
  coef2  <- par[(nrow(K) + 1):(2 * nrow(K))]
  theta <- exp(par[(2 * nrow(K)+1)])
  beta0_1 <- par[(2 * nrow(K)+2)]
  beta0_2 <- par[(2 * nrow(K)+3)]
  
  r <- -sum( 
            -log(1 + exp(beta0_1 + K%*%coef1)) + 
            (as.numeric(y != 0)) * 
              (
                lgamma(theta + y) + theta * log(theta) - 
                (theta + y)*log(theta + exp(beta0_2 + K%*%coef2)) + y*(beta0_2 + K%*%coef2) +
                (beta0_1 + K%*%coef1) - lgamma(1 + y) - lgamma(theta) - 
                log(1 - (theta / (theta + exp(beta0_2 + K%*%coef2)) )^theta)
              )
            ) + lambda1 * t(coef1)%*%K%*%coef1 + lambda2 * t(coef2)%*%K%*%coef2
  
  ## Try to avoid unrealistic values of theta
  if(theta < 1e-5 | theta > 20){
    r <- 1e13
  }
  
  return(r)
}

## The gradient for the KRHNB target function, confirmed by matching numerical derivative in R
krhnb.gr <- function(par, K, y, lambda1, lambda2){
  
  coef1  <- par[1:nrow(K)]
  coef2  <- par[(nrow(K) + 1):(2 * nrow(K))]
  theta <- exp(par[(2 * nrow(K)+1)])
  beta0_1 <- par[(2 * nrow(K)+2)]
  beta0_2 <- par[(2 * nrow(K)+3)]
  
  dc1 <- K%*%(exp(beta0_1 + K%*%coef1) / (1+exp(beta0_1 + K%*%coef1)) - as.numeric(y != 0)) +
    2*lambda1 * K%*%coef1
  
  dc2 <- K%*%( as.numeric(y != 0) * (
    ((theta + y) * exp(beta0_2 + K%*%coef2) / (theta+exp(beta0_2 + K%*%coef2))) - y +
      (theta^2 * (theta / (theta + exp(beta0_2 + K%*%coef2)))^(theta - 1) * 
         (exp(beta0_2 + K%*%coef2) / (theta+exp(beta0_2 + K%*%coef2))^2)) /
      (1 - (theta / (theta + exp(beta0_2 + K%*%coef2)))^(theta) )
  )) + 2*lambda2 * K%*%coef2
  
  dt <- theta * (-as.numeric(y != 0) %*%
                   (
                     -(theta + y) / (theta + exp(beta0_2 + K%*%coef2)) + 
                       (
                         (theta / (theta + exp(beta0_2 + K%*%coef2)))^(theta) * 
                           (log(theta / (theta + exp(beta0_2 + K%*%coef2))) + 
                              (1 - (theta/ (theta+exp(beta0_2 + K%*%coef2)))))
                       ) / 
                       (1 - (theta / (theta + exp(beta0_2 + K%*%coef2)))^(theta) ) + 
                       digamma(theta + y) + log(theta) + 1 - log(theta + exp(beta0_2 + K%*%coef2)) - digamma(theta) 
                   ))

  db0_1 <- sum(exp(beta0_1 + K%*%coef1) / (1+exp(beta0_1 + K%*%coef1)) - as.numeric(y != 0))
    
  db0_2 <- sum( as.numeric(y != 0) * (
    ((theta + y) * exp(beta0_2 + K%*%coef2) / (theta+exp(beta0_2 + K%*%coef2))) - y +
      (theta^2 * (theta / (theta + exp(beta0_2 + K%*%coef2)))^(theta - 1) * 
         (exp(beta0_2 + K%*%coef2) / (theta+exp(beta0_2 + K%*%coef2))^2)) /
      (1 - (theta / (theta + exp(beta0_2 + K%*%coef2)))^(theta) )
  ))

    return(c(dc1,dc2,dt,db0_1,db0_2))

}

## Computes the Gaussian kernel matrix from the data matrix X
gaussKernel <- function(X=NULL, sigma=NULL) {
  return( exp(-1*as.matrix(dist(X)^2)/sigma) )
}


## Optimize C and Theta in 'krhnb.fn' given the data and lambda
solveForC <- function(pars=rep(0, (2*nrow(K))+3),
                      y = NULL,
                      K = NULL,
                      lambda1 = NULL,
                      lambda2 = NULL,
                      con = list(),
                      hessian = FALSE) {
  
  chattheta <- optim(par = pars, krhnb.fn, gr=krhnb.gr, K=K, y=y, lambda1=lambda1, lambda2=lambda2,
                       method="BFGS", control=con, hessian = hessian)
  
  chat1 <- chattheta$par[1:nrow(K)]
  chat2 <- chattheta$par[(nrow(K)+1):(2*nrow(K))]
  thetahat <- exp(chattheta$par[(2*nrow(K)+1)])
  beta0_1 <- chattheta$par[(2*nrow(K)+2)]
  beta0_2 <- chattheta$par[(2*nrow(K)+3)]
  
  return(list(chat1 = chat1,
              chat2 = chat2,
              thetahat = thetahat,
              beta0_1 = beta0_1,
              beta0_2 = beta0_2,
              fullopt = chattheta))
}


## Produce predict Ys given values of c and theta and data K, using
## Parameters:
##   'c' - a nx1 vector of c values
##   'theta' - a scalar for theta
##   'K' - a kernel matrix produced by 'newKernel' function or for cross-
##         validation with dimensions [fold, -fold], or for in-sample just the
##         original kernel.
##         See the 'newKernel' function to see how you produce a Kernel for
##         prediction.
## Values:
##   'yhat' - a nx1 vector of predicted y values
predictY <- function(c1, c2, theta, beta0_1, beta0_2, K) {
  yhat <- (exp(beta0_1 + K%*%c1) * exp(beta0_2 + K%*%c2)) / 
            ((1+exp(beta0_1 + K%*%c1)) * (1 - (theta / (theta + exp(beta0_2 + K%*%c2)))^theta))
  
  return(yhat)
}



## Produce a new Kernel from new data and data that was used to fit model for
## prediction
## Parameters:
##   'X' - the original data matrix, UNSCALED
##   'newData' - the new data matrix (with same features), unscaled
##   'whichkernel' - which kernel was used on the original data
## Values:
##   'newK' - The new Kernel to be used for prediction
newKernel <- function(X, newData, whichkernel = "gaussian") {
  
  # Get values of oldData to scale newData
  Xmeans <- colMeans(X)
  Xsd <- apply(X, 2, sd)
  X <- scale(X, center = Xmeans, scale = Xsd)


  # scale test data by means and sd of training data
  newData <- scale(newData, center = Xmeans, scale = Xsd)      
  
  # predict based on new kernel matrix
  # kernel distances for test points (simply recompute all pairwise distances here because dist() is so fast )
  nn <- nrow(newData)
  
  ## Compute kernel matrix
  K <- NULL
  if(whichkernel=="gaussian"){ 
    newK <- exp(-1*as.matrix(pdist(newData, X))^2/(2*ncol(X)))
  }
  if(whichkernel=="linear"){ K <- tcrossprod(rbind(newData, X)) }
  if(whichkernel=="poly2"){ K <- (tcrossprod(rbind(newData, X))+1)^2 }
  if(whichkernel=="poly3"){ K <- (tcrossprod(rbind(newData, X))+1)^3 }
  if(whichkernel=="poly4"){ K <- (tcrossprod(rbind(newData, X))+1)^4 }
  if(is.null(K) & is.null(newK)){ stop("No valid Kernel specified") }

  if(whichkernel!="gaussian"){ 
    newK <- matrix(K[1:nn, (nn+1):(nn+nrow(X))],
                 nrow=nrow(newData),
                 byrow=FALSE)
  }
  
  return(newK)
}

## Function that splits a vector in to n chunks
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

## Function for the lambda search, returns the CV error
lambda.fn <- function(par = NULL,
                      X = NULL,
                      y = NULL,
                      theta = NULL,
                      folds = NULL,
                      chunks = NULL) {

  lambda1 <- exp(par[1])
  lambda2 <- exp(par[2])

  mse <- 0
  for(j in 1:folds){
    ## Get test set
    fold <- chunks[[j]]
    ## Get training set
    otherX <- X[-fold,,drop=FALSE]
    ## Get training K
    firstK <- gaussKernel(X = scale(otherX), sigma=ncol(otherX))
    pars <- rep(0, 2*nrow(firstK) + 3)
    chattheta <- solveForC(pars = pars, y=y[-fold], K=firstK, lambda1=lambda1, lambda2=lambda2)
    ## Get test K
    testK <- newKernel(otherX, X[fold,,drop=FALSE])
    yhat <- predictY(c1=chattheta$chat1, c2=chattheta$chat2,
                     theta=chattheta$thetahat, beta0_1 = chattheta$beta0_1,
                     beta0_2 = chattheta$beta0_2, K=testK)
    mse <- mse + sum((y[fold] - yhat)^2)
  }
  rmspe <- sqrt(mse/length(y))
  
  return(rmspe)
}

## This function returns a data frame that is the same dimensions as the
## original data. In each cell is the numerical derivative along that column
## for that observation (row).
## pars can either be a krhnb object or a solve for C object
## Must pass X.init, not transformed
effects <- function(X,
                    pars,
                    diff = 1e-5
                    #diffStages = T todo: add dummy for whether just to do overall effect
                    ){

  effectDF <- matrix(NA, ncol=ncol(X), nrow=nrow(X),
                     dimnames=list(NULL, colnames(X)))
  onesDF <- matrix(NA, ncol=ncol(X), nrow=nrow(X),
                     dimnames=list(NULL, colnames(X)))
  countDF <- matrix(NA, ncol=ncol(X), nrow=nrow(X),
                     dimnames=list(NULL, colnames(X)))
  
  X <- scale(X)
  
  for(i in 1:ncol(X)){
    
    testX1 <- X
    testX2 <- X
    binaryindicator <- FALSE
    if (length(unique(X[, colnames(X)[i]])) == 2){
      Xquant <- unique(X[, colnames(X)[i]]) # if dichotomous, just choose the two values
      testX1[, colnames(X)[i]] <- Xquant[1]
      testX2[, colnames(X)[i]] <- Xquant[2]
      binaryindicator <- TRUE
    } else{
      Xquant <- cbind(X[, colnames(X)[i]], X[, colnames(X)[i]] + diff)
      testX1[, colnames(X)[i]] <- Xquant[,1]
      testX2[, colnames(X)[i]] <- Xquant[,2]
    }
    
    nK1 <- newKernel(X, testX1) # todo: add whichkernel
    nK2 <- newKernel(X, testX2)
    
    yhat1 <- predictY(c1=pars$chat1, c2=pars$chat2, theta=pars$thetahat,
                      beta0_1 = pars$beta0_1, beta0_2 = pars$beta0_2, K=nK1)
    yhat2 <- predictY(c1=pars$chat1, c2=pars$chat2, theta=pars$thetahat,
                      beta0_1 = pars$beta0_1, beta0_2 = pars$beta0_2, K=nK2)
    
    ydif <- yhat2 - yhat1
    if (!binaryindicator) ydif <- ydif/diff
    effectDF[, i] <- ydif
    
    pr1.1 <- predict1(c1=pars$chat1, c2=pars$chat2, theta=pars$thetahat,
                      beta0_1 = pars$beta0_1, beta0_2 = pars$beta0_2, K=nK1)
    pr1.2 <- predict1(c1=pars$chat1, c2=pars$chat2, theta=pars$thetahat,
                      beta0_1 = pars$beta0_1, beta0_2 = pars$beta0_2, K=nK2)
    
    pr1dif <- pr1.2 - pr1.1
    if (!binaryindicator) pr1dif <- pr1dif/diff
    
    onesDF[, i] <- pr1dif
    
    mu1 <- predictmu(c2=pars$chat2, theta = pars$thetahat, beta0_2 = pars$beta0_2, K=nK1)
    mu2 <- predictmu(c2=pars$chat2, theta = pars$thetahat, beta0_2 = pars$beta0_2, K=nK2)
    
    mudif <- mu2 - mu1
    if (!binaryindicator) mudif <- mudif/diff
    countDF[, i] <- mudif
  }
  return(list(effectDF=effectDF,
              onesDF=onesDF,
              countDF=countDF))
}

## To plot the histograms of the krhnb object, under development
plot.krhnb <- function(krhnb.obj, type = "effect", type2 = "both", varNames = NULL, zeroline = TRUE){

  if(type == "median"){
    mat <- krhnb.obj$bootmfx.median
  } else if(type == "mean"){
    if (type2 == "both") {
      mat <- krhnb.obj$bootmfx
      meaneff <- krhnb.obj$meaneffect
    } else if (type2 == "ones") {
      mat <- krhnb.obj$bootonesmfx
      meaneff <- apply(krhnb.obj$onesmat, 2, mean)
    } else if (type2 == "count") {
      mat <- krhnb.obj$bootcountmfx
      meaneff <- apply(krhnb.obj$countmat, 2, mean)
    }

    quantiles <- apply(mat, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
    plot.df <- data.frame(variable = colnames(mat),
                          mean.est = meaneff,
                          p2.5 = quantiles[1, ],
                          p97.5 = quantiles[2, ])
    plot.df$variable <- reorder(plot.df$variable, meaneff)

    p <- ggplot(plot.df, aes(y=mean.est, x=as.factor(variable), ymin = p2.5, ymax = p97.5))
    if(zeroline) {
      p <- p + geom_hline(yintercept = 0, color = "red", size = 1.1)
    }
    p <- p + geom_pointrange(size = 0.8) +
      coord_flip() +  
      theme_bw()

    
  } else if(type == "effect") {
    
    mat <- k.out$effectmat
    medians <- apply(mat, 2, median)
    plot.df <- melt(mat[, names(sort(medians))])
    
    p <- ggplot(plot.df, aes(y=value, x=Var2)) +
      geom_violin(fill = "lightblue") +
      coord_flip() +  
      theme_bw()
  }
  
  p
}

## Give CIs on average pwmfx
summary.krhnb <- function(krhnb.obj, type = "median"){
  
  if(type == "median"){
    mat <- krhnb.obj$bootmfx.median
  } else if(type == "mean"){
    mat <- krhnb.obj$bootmfx.mean
  } else{
    stop('Only supports type "median" or "mean" right now.')
  }
  
  sapply(colnames(mat),
         function(x) quantile(mat[, x],
                              probs = c(0.025, 0.5, 0.975)))
}

## Predict 1s from logit, first component
predict1 <- function(c1, c2, beta0_1, beta0_2, theta, K) {
  
  pry1 <- 1/ (1+exp(-(beta0_1 + K%*%c1)))
  return(pry1)
}

## Predict mean of truncated negative binomial, second component
predictmu <- function(c2, theta, beta0_2, K) {
  
  mu <- exp(beta0_2 + K%*%c2) / 
    pnbinom(0, size = theta, mu = exp(beta0_2 + K%*%c2), lower.tail = FALSE)  
  # strictly speaking, not calculating mu (=E[nbinom]), but E[truncated nbinom] = E[nbinom]/pnbinom(>0) ; good
  return(mu)
}