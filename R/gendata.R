gendata.bayes <- function(n0,n1,m0,m1,p,alpha)
{
	theta.c0 = numeric(p)
	theta.c1 = numeric(p)
	
	P = runif(p,0,1)
	
	for(i in 1:p)
	{	theta.c0[i]=rbeta(1,alpha*P[i],alpha*(1-P[i]))
		theta.c1[i]=rbeta(1,alpha*P[i],alpha*(1-P[i]))
	}		
	train = matrix(runif(n0*p+n1*p),n0+n1,p)
	test = matrix(runif(m0*p+m1*p),m0+m1,p)
			
	train = 1* (train < rbind(matrix(rep(theta.c0,n0),n0,p,TRUE),
            	        	  matrix(rep(theta.c1,n1),n1,p,TRUE)))
	test = 1* (test < rbind(matrix(rep(theta.c0,m0),m0,p,TRUE),
            	        	matrix(rep(theta.c1,m1),m1,p,TRUE)))
	train <- cbind(c(rep(0,n0),rep(1,n1)),train)
	test <- cbind(c(rep(0,m0),rep(1,m1)),test)
	list(train=train,test=test)
}
