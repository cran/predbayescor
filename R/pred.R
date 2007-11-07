################################################################################
#this function makes prediction on test data 
#  test --- the test data, the first column is response, which could be NA 
#           if not available
#  train --- the training data set, the first column is the response
#  is.binary.features --- the indicator whether the features are binary
#  subset.sel --- the index in train data used to select features
#  k --- the number of features retained
#  theta0 --- the left bound of the prior of theta, 
#  the right bound is required to be 1-theta0 
#  no.theta --- the number of theta used in simpson rule
#  no.alpha --- the number of alpha's, which is determined with the quantiles of
#               inverse gamma distribution defined with
#  alpha.shape --- the shape parameter and 
#  alpha.rate  --- the rate parameter 
#  correct --- the indicator whether the correction method shall be applied
#  no.theta.adj  --- the number of thetas used in adjustment factor

Predict.bayes <- function(test,train,is.binary.features=FALSE,k,
                          subset.sel=1:nrow(train),
                          theta0=0,no.theta=20,
                          alpha.shape=0.5,alpha.rate=5,no.alpha=5,
                          correct=TRUE,no.theta.adj=20)
{
    #convert the real-values into binary
    if(!is.binary.features)
    {  thresholds <- find.thresholds(rbind(train[,-1],test[,-1]))
       train[,-1] <- real2bin(train[,-1],thresholds)
    }
    p <- ncol(train) - 1
    no.test <- nrow(test)
    no.train <- nrow(train)
    no.train0 <- sum(train[,1]==0)
    no.train1 <- no.train - no.train0
    theta.range <- c(theta0,1-theta0)
    #selecte features with absoluate correlation
    if(k < p)
    {   info.sel <- selfth.abscor(k,x=train[subset.sel,-1],
                     y=train[subset.sel,1])
        fth.sel <- info.sel$fth+1
        gamma.sel <- info.sel$abscors[k]
    
    }
    else 
    {  fth.sel <- 1:p+1
       gamma.sel <- 0
    } 
    #converting real-values of selected features in test into binary
    if(!is.binary.features) 
       test[,fth.sel] <- real2bin(test[,fth.sel],thresholds[fth.sel-1])
    #discard the unselected features to save memory   
    train <- train[,c(1,fth.sel),drop=FALSE]
    test <- test[,c(1,fth.sel),drop=FALSE] 
    #the number of ones for selected features
    I1.s <- apply(train[train[,1]==1,-1,drop=FALSE],2,sum)
    I0.s <- apply(train[train[,1]==0,-1,drop=FALSE],2,sum)
    #make the vector of all possible alphas
    alpha.bank <- gen.alpha.IGamma(no.alpha,alpha.shape,alpha.rate) 
    log.prior <- dgamma(1/alpha.bank,alpha.shape,alpha.rate,log=TRUE) - 
                 2*log(alpha.bank)     		   
    #if correction is required, calcalte the adjustment factor over alpha.bank
    if(correct & k < p)
    {    #construct select region which is a list of two elements, 
         #I0:a vector,the upper bounds of I0 for each value of I1
         #b: a number, the smallest value of I1
	  #number of class 0 used to select features
         n0 <- sum(train[subset.sel,1]== 0) 
         n1 <- length(subset.sel) - n0 
         y.bar <- n1/(n0+n1)
         b <- ceiling(1/(1/length(subset.sel)+(1-y.bar)/(n1*gamma.sel^2)))	
         S.gamma <- .C("select_bd",as.integer(n0),as.integer(n1),gamma.sel,
                  b=as.integer(b),I0=as.integer(rep(0,n1-b+1)),
		  PACKAGE="predbayescor")[4:5]
         adj.factor <- .C("prob_lowcor",as.integer(no.theta.adj),
		           as.integer(n0),as.integer(n1),as.integer(no.alpha),
		           length(S.gamma$I0),as.integer(S.gamma$I0),
			   as.integer(S.gamma$b),as.integer(max(S.gamma$I0)),
			   alpha.bank,theta.range,
			   probs.lowcor=rep(0,no.alpha),
			   PACKAGE="predbayescor")$probs.lowcor	
    }
    else adj.factor <- rep(1,no.alpha)
    #make the log.prior over alpha.bank which has been adjusted with
    #adjustment factor
    log.adj.factor <- log(adj.factor) * (p-k)        
            
    #call C function to do the prediction
    out.c <- .C("predict_int_c",as.integer(k),as.integer(no.test),
                   as.integer(no.train0),as.integer(no.train1),
		   as.integer(no.alpha),as.integer(no.theta),
		   as.integer(t(test[,-1])),theta.range,
		   as.integer(I1.s),as.integer(I0.s),log.adj.factor, 
		   alpha.bank, prob.pred=rep(0,no.test),
                   as.integer(1),log.margin=rep(0,no.alpha),
		   PACKAGE="predbayescor")
    prob.pred <- out.c$prob.pred
    log.post.alpha <- out.c$log.margin + log.prior
    alpha.prior.adj.post <- cbind(alpha.bank,log.prior,
                                  adj.factor,log.post.alpha)	   
    #summarize the prediction
    #predict
    states.pred <- 1*(prob.pred>0.5) 
    wrong <- 1*(states.pred != test[,1])
    error.rate <- mean(wrong) 
    # calculate average minus log probs
    aml <- -mean(dbinom(test[,1],1,prob.pred,log=TRUE))
    # tabulate the predictive probs
    summary.pred <- present(prob.pred,test[,1]) 
    # calculate mean square error
    mse <- mean((prob.pred-test[,1])^2)   
    #output
    list(aml=aml,error.rate=error.rate,mse=mse,
         summary.pred=summary.pred,
         prediction=cbind(true=test[,1],pred=states.pred,
	                  prob.pred=prob.pred,wrong),
         alpha.prior.adj.post=alpha.prior.adj.post,
	  features.selected = c(feature=fth.sel,cutoff=gamma.sel))  
}
################################################################################
#calculating the absolute correlation between x and y
abs_cor <- function(x,y)
{
    if(sum(abs(x-mean(x)))==0 || sum(abs(y-mean(y)))==0 )
	return(0)
    else
	return(abs(cor(y,x)))
}
#this function selects features from x by the absolate correlation with y
selfth.abscor <- function(k,x,y)
{  abscors <- apply(x,2,abs_cor,y)
   info.sort <- sort(abscors,decreasing=TRUE,index.return=TRUE)
   list(fth=info.sort$ix[1:k], abscors=info.sort$x[1:k])
}
################################################################################
#this function generates no.alpha alpha's from Inverse Gamma distribution with
#shape and rate parameters, by taking the quantiles corresponding to the 
#the probabilities spaced equally from w to 1-w, where w = 1/no.alpha/2
gen.alpha.IGamma <- function(no.alpha,shape,rate)
{    w <- 1/no.alpha/2
     1/qgamma(seq(1 - w, w, length = no.alpha),shape,rate)

}
################################################################################

#this function converts the data, a matrix, into binary data by 
#thresholding at "thresholds"
real2bin <- function(data, thresholds=find.thresholds(data))
{  (data > matrix(thresholds,nrow(data),ncol(data),byrow=TRUE))*1   
}
#this function finds the thresholds of features, either the median for 
#real values or the mean of binary values, for the usage in function
#"real2bin"
find.thresholds <- function(data)
{   thresholds <- rep(0,ncol(data))
    for(i in 1:ncol(data))
    {   if(length(unique(data[,i]))==2)
            thresholds[i] <- mean(data[,i])
        else thresholds[i] <- median(data[,i]) 
    }
    thresholds   
}
