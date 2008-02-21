###########################################################################
cv.bayes <- function(data,is.binary.features=FALSE,no.folds=10,k,
                     theta0=0,no.theta=20,
                     alpha.shape=0.5,alpha.rate=5,no.alpha=5,
                     correct=TRUE,no.theta.adj=20)
{  no.cases <- nrow(data)
   m <- floor(no.cases/no.folds)
   rm <- no.cases - m*no.folds
   prob.pred <- rep(0,no.cases)
   features <- data.frame(feature=matrix(0,no.folds,k),cutoff=rep(0,no.folds))
   for( n in 1:no.folds ){
      if(n <= rm)
         test.range <- n + seq(0,m,by=1)*no.folds
      else test.range <- n + seq(0,m-1,by=1)*no.folds
      result.pred <-
         predict_bayes(data[test.range,,drop=FALSE],
	   data[-test.range,,drop=FALSE],
	   is.binary.features,k,,theta0,no.theta,alpha.shape,
	   alpha.rate,no.alpha,correct,no.theta.adj)
      prob.pred[test.range] <- result.pred$prediction[,"prob.pred"]
      features[n,] <- result.pred$features
   }
   #summarize the prediction
   #predict
   states.pred <- 1*(prob.pred>0.5)
   wrong <- 1*(states.pred != data[,1])
   error.rate <- mean(wrong)
   # calculate average minus log probs
   aml <- -mean(dbinom(data[,1],1,prob.pred,log=TRUE))
   # calculate mean square error
   mse <- mean((prob.pred-data[,1])^2)
   # tabulate the predictive probs
   summary.pred <- present(prob.pred,data[,1])

   #output
   list(aml=aml,error.rate=error.rate,mse=mse,
         summary.pred=summary.pred,
         prediction=cbind(true=data[,1],pred=states.pred,wrong=wrong,
                          prob.pred=prob.pred),
	 features=features)

}
