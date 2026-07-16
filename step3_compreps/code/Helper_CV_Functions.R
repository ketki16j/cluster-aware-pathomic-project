### Author: Jeremy Rubin
## Date: 1/30/26
## Helper functions for cross-validation related things

# Make folds for internal five-fold cross-validation based on the training data
CV.make.folds = function(nTrain)
{
  # Compute sample sizes
  sampleSizeF <- floor(1/5 * nTrain)
  
  # Create the randomly-sampled indices for the dataframe. Use setdiff() to avoid 
  # overlapping subsets of indices
  indicesF1 <- sort(sample(1:nTrain,size=sampleSizeF))
  indicesNotF <- setdiff(1:nTrain,indicesF1)
  
  indicesF2 <- sort(sample(indicesNotF, size=sampleSizeF))
  indicesNotF <- setdiff(1:nTrain,c(indicesF1,indicesF2))
  
  indicesF3 <- sort(sample(indicesNotF, size=sampleSizeF))
  indicesNotF <- setdiff(1:nTrain,c(indicesF1,indicesF2,indicesF3))
  
  indicesF4 <- sort(sample(indicesNotF, size=sampleSizeF))
  indicesF5 <- setdiff(indicesNotF,indicesF4)
  
  return(list(indicesF1,indicesF2,indicesF3,indicesF4,indicesF5))
}

# Uses cross-validation to get internally-validated MSE/C-statistic for 
# given algorithm and continuous or binary outcome
CV.avg.all.folds = function(X.train,Y.train,algorithm,outcome,lambda){
  
  nfolds <- 5
  
  if(algorithm!="KDPI"){folds <- CV.make.folds(length(Y.train))}
  else{folds <- CV.make.folds(length(X.train))}
  
  if(outcome=="binary")
  {
    ### Code for binary outcomes to ensure that each fold has both cases and controls 
    fold1.outcomes <- Y.train[folds[[1]]]
    fold2.outcomes <- Y.train[folds[[2]]]
    fold3.outcomes <- Y.train[folds[[3]]]
    fold4.outcomes <- Y.train[folds[[4]]]
    fold5.outcomes <- Y.train[folds[[5]]]
    
    while(length(unique(fold1.outcomes)) == 1 | 
          length(unique(fold2.outcomes)) == 1 | 
          length(unique(fold3.outcomes)) == 1 | 
          length(unique(fold4.outcomes)) == 1 | 
          length(unique(fold5.outcomes)) == 1)
    {
      print(paste("Redoing folds for cross-validation because not all folds have both a case and control"))
      if(algorithm!="KDPI"){folds <- CV.make.folds(length(Y.train))}
      else{folds <- CV.make.folds(length(X.train))}
      
      fold1.outcomes <- Y.train[folds[[1]]]
      fold2.outcomes <- Y.train[folds[[2]]]
      fold3.outcomes <- Y.train[folds[[3]]]
      fold4.outcomes <- Y.train[folds[[4]]]
      fold5.outcomes <- Y.train[folds[[5]]]
    }  
  }
  
  CV.fold.matrix <- diag(nfolds)
  mse.vec <- rep(0,nfolds)
  
  ## Iterate through different train/test fold setups
  for(i in 1:length(folds))
  {
    fold.setup <- CV.fold.matrix[i,]
    train.fold <- which(fold.setup==0)
    test.fold <- which(fold.setup==1)
    cv.train.ind <- sort(c(unlist(folds[train.fold])))
    cv.test.ind <- sort(c(unlist(folds[test.fold])))
    
    ## KDPI uses just a single vector/predictor whereas other machine learning models
    ## use multiple predictors in the model fitting process
    if(algorithm!="KDPI")
    {
      if(is.vector(X.train) && !is.matrix(X.train) && !is.data.frame(X.train))
      {
        mse.vec[i] <- CV.mse.one.setup(X.train[cv.train.ind],
                                       Y.train[cv.train.ind],
                                       X.train[cv.test.ind],
                                       Y.train[cv.test.ind],
                                       algorithm,
                                       outcome,
                                       lambda)
      }
      else{mse.vec[i] <- CV.mse.one.setup(X.train[cv.train.ind,],
                                          Y.train[cv.train.ind],
                                          X.train[cv.test.ind,],
                                          Y.train[cv.test.ind],
                                          algorithm,
                                          outcome,
                                          lambda)}
        
    } else{
      
      mse.vec[i] <- CV.mse.one.setup(X.train[cv.train.ind],
                                     Y.train[cv.train.ind],
                                     X.train[cv.test.ind],
                                     Y.train[cv.test.ind],
                                     algorithm,
                                     outcome,
                                     lambda)
      }
    
  }
  
  return(mean(mse.vec))
  
}

# Compute MSE/C-statistic for given training/testing fold setup
# algorithm refers to whether you're looking at the lasso, ridge regression, 
# elastic net, or random forest
CV.mse.one.setup = function(X.train,Y.train,X.test,Y.test,algorithm,outcome,lambda)
{
  # eGFR outcome prediction
  if(outcome=="continuous")
  {
      if(algorithm=="lasso")
      {
        mod <- glmnet(x=data.matrix(X.train),y=Y.train,alpha=1)
        predictions <- predict(mod, newx=data.matrix(X.test),s=lambda)
      } 
      
      if(algorithm=="ridge")
      {
        mod <- glmnet(x=data.matrix(X.train),y=Y.train,alpha=0)
        predictions <- predict(mod, newx=data.matrix(X.test),s=lambda)
      }
      
      if(algorithm=="enet")
      {
        mod <- glmnet(x=data.matrix(X.train),y=Y.train,alpha=0.5)
        predictions <- predict(mod, newx=data.matrix(X.test),s=lambda)
      }
      
     if(algorithm=="KDPI"){
       
       # Create proper dataframes for both training and testing
       train.df <- data.frame(X_train = X.train, Y_train = Y.train)
       test.df <- data.frame(X_train = X.test)  # Use the same variable name as training
       
       # Fit the model
       glm.mod <- glm(Y_train ~ X_train, data=train.df, family=gaussian(link="identity"))
       
       # Make predictions - note we're using the correct variable name in test.df
       predictions <- predict(glm.mod, newdata=test.df, type="response")
    }
    
    if(algorithm=="rf") {
      
      if(is.vector(X.train) && !is.matrix(X.train) && !is.data.frame(X.train))
      {
        train.df <- cbind(0,X.train)
        colnames(train.df) <- c("zero","X.train")
        
        mod <- randomForest(
          x = train.df,  
          y = Y.train,                                              
          ntree = 500,
          nodesize=lambda,
          keep.forest=TRUE,
          keep.inbag=TRUE
        )  
        
        test.df <- cbind(0,X.test)
        
        # Use the same variable name as training
        colnames(test.df) <- c("zero","X.train")
        
        predictions <- predict(mod,newdata=test.df)
        
      } else{
        mod <- randomForest(
          x = X.train,  
          y = Y.train,                                              
          ntree = 500,
          nodesize=lambda,
          keep.forest=TRUE,
          keep.inbag=TRUE
        )
        
        predictions <- predict(mod,newdata=X.test)
      }
    }
    
    MSE <- mean((Y.test - predictions)^2)
    return(MSE)
      
  }
  
  # DGF prediction
  if(outcome=="binary")
  {
      if(algorithm=="lasso")
      {
        mod <- glmnet(x=data.matrix(X.train),y=Y.train,alpha=1,family="binomial",type.measure="auc")
        C.statistic <- assess.glmnet(mod,newx=as.matrix(X.test),newy=Y.test,family="binomial",s=lambda)[[3]][1]
      } 
      
      if(algorithm=="ridge")
      {
        mod <- glmnet(x=data.matrix(X.train),y=Y.train,alpha=0,family="binomial",type.measure="auc")
        C.statistic <- assess.glmnet(mod,newx=as.matrix(X.test),newy=Y.test,family="binomial",s=lambda)[[3]][1]
        
      }
      
      if(algorithm=="enet")
      {
        mod <- glmnet(x=data.matrix(X.train),y=Y.train,alpha=0.5,family="binomial",type.measure="auc")
        C.statistic <- assess.glmnet(mod,newx=as.matrix(X.test),newy=Y.test,family="binomial",s=lambda)[[3]][1]
        
      }
      
      if(algorithm=="rf")
      {
        if(is.vector(X.train) && !is.matrix(X.train) && !is.data.frame(X.train))
        {
          train.df <- cbind(0,X.train)
          colnames(train.df) <- c("zero","X.train")
          
          mod <- randomForest(
            x = train.df,  
            y = Y.train,                                              
            ntree = 500,
            nodesize=lambda,
            keep.forest=TRUE,
            keep.inbag=TRUE
          )  
          
          test.df <- cbind(0,X.test)
          colnames(test.df) <- c("zero","X.train")
          # Use the same variable name as training
          
          predictions <- predict(mod,newdata=test.df,type="prob")
          
        } else{
          mod <- randomForest(
            x = X.train,  
            y = Y.train,                                              
            ntree = 500,
            nodesize=lambda,
            keep.forest=TRUE,
            keep.inbag=TRUE
          )
          
          predictions <- predict(mod,newdata=X.test,type="prob")
        }
        
        
        C.statistic <- auc(Y.test,predictions[,2])
        
      }
    
      
    if(algorithm=="KDPI"){
      
      # Can compute the C-statistic directly using the KDPI score, no model 
      # fitting needed
      C.statistic <- auc(Y.test,X.test)
    }
      
      return(C.statistic)
      
    }
  }

# Return MSE for given lambda
glmnet.opt.lambda <- function(lambda, glmnet.mod, newx, newy) {
  
  # Get predictions for the specific lambda value
  preds <- predict(glmnet.mod, newx = newx, s = lambda)
  
  # Calculate MSE
  mse <- mean((newy - preds)^2)
  
  return(mse)
}

# Return C-statistic for given lambda
glmnet.opt.lambda.binary <- function(lambda, glmnet.mod, newx, newy)
{
  C.statistic <- assess.glmnet(glmnet.mod,newx=as.matrix(newx),
                              newy=newy,
                              family="binomial",
                              s=lambda)[[3]][1]
  return(C.statistic)
}
