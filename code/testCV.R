compareKrhnb <- function(y, X, folds, hnb = TRUE, ...) {
  
  ## Container DFs
  totalPredDF <- NULL
  totalMseDF <- NULL
  
  foldids <- 1:length(folds)
  
  for(i in foldids){
    training <- setdiff(1:nrow(X), folds[[i]])
    ## Create the training data using these indices
    trainingX <- X[training, ]
    trainingy <- y[training]
    
    ## Create the test data using the ith set of indices
    testX <- X[folds[[i]], ]
    testy <- y[folds[[i]]]
    
    out <- predictFits(testy, trainingy,
                       testX, trainingX,
                       i, hnb = hnb, ...)
    
    totalPredDF <- rbind(totalPredDF, out$predictions)
    totalMseDF <- rbind(totalMseDF, out$mse)
  }

  return(list(totalPredDF = totalPredDF,
              totalMseDF = totalMseDF))
}

## Do one iteration of the cross validation
predictFits <- function(testy, trainingy,
                        testX, trainingX,
                        i, hnb, ...){
  ## Create the lists of training indices
  
  ## Generate dichotomous outcomes for the two-stage model
  testybin <- as.numeric(testy!=0)
  trainingybin <- as.numeric(trainingy!=0)
  
  ## Generate data frames for some of the methods that are a little particular
  trainingDFbin <- cbind(trainingybin, trainingX)
  trainingDF <- cbind(trainingy, trainingX)
  testDFbin <- cbind(testybin, testX)
  testDF <- cbind(testy, testX)
  
  ## Fitting KRLS
  k <- krls(X=trainingX, y=trainingy)
  yhat.k <- predict(k, testX)$fit
  #print(yhat.k)
  yhat.trunc0.k <- ifelse(yhat.k < 0, 0, yhat.k)
  #print(yhat.k)
  mse.k <- mean( (testy - yhat.k)^2 )
  mae.k <- mean( abs(testy - yhat.k) )
  
  mse.trunck <- mean( (testy - yhat.trunc0.k)^2 )
  mae.trunck <- mean( abs(testy - yhat.trunc0.k) )
  
  ## Fitting KRHNB
  k.nb <- krhnb(X=trainingX, y=trainingy, bootstrap = F, ...)
  
  testKern <- newKernel(X=trainingX, newData=testX)
  yhat.khnb <- predictY(c1 = k.nb$chat1,
                       c2 = k.nb$chat2,
                       theta = k.nb$thetahat,
                       beta0_1=k.nb$beta0_1,
                       beta0_2=k.nb$beta0_2,
                       K = testKern)
  mse.knb <- mean( (testy - as.numeric(yhat.khnb))^2 )
  mae.knb <- mean( abs(testy - as.numeric(yhat.khnb)) )
  
  ## Hurdle
  if(hnb){
  hnb <- hurdle(trainingy ~ .,
                data = trainingDF,
                dist = "negbin",
                zero.dist = "binomial",
                link = "logit")
                
  yhat.hnb <- predict(hnb, testDF[, 2:ncol(testDF)], type = "response")
  mse.hnb <- mean( (testy - as.numeric(yhat.hnb))^2 )
  mae.hnb <- mean( abs(testy - as.numeric(yhat.hnb)) )
  } else {
    yhat.hnb <- NA
    mse.hnb <- NA
    mae.hnb <- NA
  }
  
  ## OLS
  lmout <- lm(trainingy ~ .,
                  data = trainingDF)
  yhat.lm <- predict(lmout, testDF[, 2:ncol(testDF)])
  mse.lm <- mean( (testy - as.numeric(yhat.lm))^2 )
  mae.lm <- mean( abs(testy - as.numeric(yhat.lm)) )
  
  ## Random Forest
  rf <- randomForest(x=trainingX, y=trainingy)
  # oos predictions for:
  # randomforest
  yhat.rf <- predict(rf, testX)
  
  mse.rf <- mean( (testy - yhat.rf)^2 )
  mae.rf <- mean( abs(testy - yhat.rf) )
  
  return(list(predictions=data.frame(testy = testy,
                                     yhat.k = yhat.k,
                                     yhat.trunc0.k = yhat.trunc0.k,
                                     yhat.knb = yhat.khnb,
                                     yhat.lm = yhat.lm,
                                     yhat.rf = yhat.rf,
                                     fold = i),
              mse=data.frame(mse.k = mse.k,
                             mae.k = mae.k,
                             mse.trunck = mse.trunck,
                             mae.trunck = mae.trunck,
                             mse.knb = mse.knb,
                             mae.knb = mae.knb,
                             mse.lm = mse.lm,
                             mae.lm = mae.lm,
                             mse.rf = mse.rf,
                             mae.rf = mae.rf,
                             fold = i)))
}
