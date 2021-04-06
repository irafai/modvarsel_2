function_par<-function(X,Y,train,method,nperm=nperm, cutoff = cutoff,nbsel = nbsel,ntree=ntree){
  
  Xtrain <- X[train,]
  Ytrain <- Y[train]
  
  imp <- varimportance(Xtrain, Ytrain, method = method,nperm=nperm,ntree=ntree)
  selvar <- select(imp, cutoff = cutoff,nbsel = nbsel)
  
  Xtest <- X[-train,]
  Ytest <- Y[-train]
  
  
  Xtest_sel <- Xtest[,selvar$indices, drop=FALSE]
  Xtrain_sel <- Xtrain[,selvar$indices, drop=FALSE]
  varsel <- selvar$var
  
  #========
  # RegLin
  #========
  if (method=="linreg"){
    model <- stats::lm(Ytrain~.,data = data.frame(Xtrain_sel))
    Ypred <- stats::predict(model,newdata = data.frame(Xtest_sel))
    
    mse <- mean((Ytest-Ypred)^2)
    
    #with all available variables
    model <- stats::lm(Ytrain~.,data = data.frame(Xtrain))
    Ypred_c <- stats::predict(model,newdata = data.frame(Xtest))
    mse_c <- mean((Ytest-Ypred_c)^2)
  }    
  #========
  # SIR
  #========
  else if (method=="sir"){  
    beta <- edrGraphicalTools::edr(Ytrain, Xtrain_sel,
                                   H = 10, K = 1, method = "SIR-I")$matEDR[, 1,
                                                                           drop = FALSE]
    indice <- Xtrain_sel%*%beta
    hopt <- cv_bandwidth(indice, Ytrain,graph.CV=FALSE)$hopt
    indicetest <- Xtest_sel%*%beta
    mat <- cbind(Ytest,indicetest)
    matord <- mat[order(mat[,2]),]
    Ytest_ord <- matord[,1]
    indicetest_ord <- matord[,2]
    Ypred <- stats::ksmooth(indice, Ytrain, kernel = "normal",
                            bandwidth = hopt,
                            x.points = indicetest_ord)$y
    
    mse <- mean((Ytest_ord-Ypred)^2)
    varsel <- selvar$var
    
    #with all available variables
    beta <- edrGraphicalTools::edr(Ytrain, Xtrain,
                                   H = 10, K = 1, method = "SIR-I")$matEDR[, 1,
                                                                           drop = FALSE]
    indice <- Xtrain%*%beta
    hopt <- cv_bandwidth(indice, Ytrain,graph.CV=FALSE)$hopt
    indicetest <- Xtest%*%beta
    mat <- cbind(Ytest,indicetest)
    matord <- mat[order(mat[,2]),]
    Ytest_ord <- matord[,1]
    indicetest_ord <- matord[,2]
    Ypred <- stats::ksmooth(indice, Ytrain, kernel = "normal",
                            bandwidth = hopt,
                            x.points = indicetest_ord)$y
    
    mse_c <- mean((Ytest_ord-Ypred)^2)
  }
  #================
  # Random Forests
  #================
  else if (method=='rf'){
    #with variables selection
    model <- randomForest::randomForest(x = Xtrain_sel,y = Ytrain, ntree = ntree)
    Ypred2 <-stats::predict(model,newdata = Xtest_sel,type = "response")
    
    mse <- mean((Ytest-Ypred)^2)
    varsel <- selvar$var
    
    #with all available variables
    model <- randomForest::randomForest(x = Xtrain,
                                        y = Ytrain, ntree = ntree)
    Ypred <- stats::predict(model,newdata = Xtest,
                            type = "response")
    
    mse_c[i] <- mean((Ytest-Ypred)^2)
  }
  
  #========
  # PCR
  #========
  else if (method=="pcr"){
    #with variables selection
    model <- pls::pcr(Ytrain~., data = data.frame(Xtrain_sel),
                      validation = "CV", scale=FALSE)
    # rmsep <- pls::RMSEP(model, intercept = FALSE)$val["CV",,]
    rmsep <- sqrt(model$validation$PRESS/n_train)
    ncomp <- find.cpt(rmsep)
    Ypred <- as.vector(stats::predict(model, data.frame(Xtest_sel),
                                      ncomp=ncomp))
    mse <- mean((Ytest - Ypred)^2)
    varsel <- selvar$var
    
    #with all available variables
    model <- pls::pcr(Ytrain~., data = as.data.frame(Xtrain),
                      validation="CV", scale=FALSE)
    # rmsep <- pls::RMSEP(model, intercept = FALSE)$val["CV",,]
    rmsep <- sqrt(model$validation$PRESS/n_train)
    ncomp <- find.cpt(rmsep)
    Ypred <-  as.vector(stats::predict(model, data.frame(Xtest),
                                       ncomp=ncomp))
    
    mse <- mean((Ytest - Ypred)^2)
  }
  #========
  # PLSR
  #========
  else if (method == "plsr" ){
    #with variables selection
    model <- pls::plsr(Ytrain~., data = data.frame(Xtrain_sel),
                       validation = "CV", scale=TRUE)
    # rmsep <- pls::RMSEP(model, intercept = FALSE)$val["CV",,]
    rmsep <- sqrt(model$validation$PRESS/n_train)
    ncomp <- find.cpt(rmsep)
    Ypred <- as.vector(stats::predict(model, data.frame(Xtest_sel),ncomp=ncomp))
    
    mse <- mean((Ytest - Ypred)^2)
    varsel <- selvar$var
    
    #with all available variables
    model <- pls::plsr(Ytrain~., data = as.data.frame(Xtrain),
                       validation="CV", scale=TRUE)
    # rmsep <- pls::RMSEP(model, intercept = FALSE)$val["CV",,]
    rmsep <- sqrt(model$validation$PRESS/n_train)
    ncomp <- find.cpt(rmsep)
    Ypred <-  as.vector(stats::predict(model, data.frame(Xtest),
                                       ncomp=ncomp))
    
    mse_c<- mean((Ytest - Ypred)^2)
  }
  #========
  # RIDGE
  #========
  else if (method == "ridge"){
    model <- glmnet::cv.glmnet(Xtrain_sel, Ytrain,
                               family = "gaussian", alpha=0, standardize = FALSE,
                               nfolds = 10, grouped=FALSE)
    Ypred <- stats::predict(model, Xtest_sel, s = "lambda.min")
    mse <- mean((Ytest-Ypred)^2)
    varsel <- selvar$var
    
    #with all available variables
    model <- glmnet::cv.glmnet(Xtrain, Ytrain,
                               family = "gaussian", alpha=0, standardize = FALSE,
                               nfolds = 10, grouped=FALSE)
    Ypred <- stats::predict(model, Xtest, s = "lambda.min")
    mse_c<- mean((Ytest-Ypred)^2)
  }
  list(mse,mse_c,varsel)
}
