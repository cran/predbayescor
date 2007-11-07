\name{Predict.bayes}

\alias{Predict.bayes}
\alias{cv.bayes}

\title{Classification rule based on Bayesian 
naive Bayes models with feature selection bias corrected
}

\description{\code{Predict.bayes} predicts the binary response based on  high
dimemsional binary features with Bayesian naive Bayes models. It also accepts real
values but they will be converted into  binary by thresholding at the medians
estimated from the data.  A smaller number of features can be selected based on the
correlations  with the response. The bias due to the selection procedure can be
corrected. \code{cv.bayes} is the short-cut function for cross-validation with 
\code{Predict.bayes}.  }

\usage{
Predict.bayes(test,train,is.binary.features=FALSE,k,
              subset.sel=1:nrow(train),
              theta0=0,no.theta=20,
              alpha.shape=0.5,alpha.rate=5,no.alpha=5,
              correct=TRUE,no.theta.adj=20)

cv.bayes(data,is.binary.features=FALSE,no.folds=10,k,
          theta0=0,no.theta=20,
          alpha.shape=0.5,alpha.rate=5,no.alpha=5,
          correct=TRUE,no.theta.adj=20) 
}
\arguments{

\item{test}{the test data, a matrix, i.e. the data for which we want to predict
            the responses. The row stands for the cases. 
           The first column is the binary response, 
	   which could be NA if they are missing.}

\item{train}{the training data, of the same format as test}

\item{data}{the data used in cross-validation, of the same format as test}

\item{no.folds}{the number of blocks the data is divided into in
cross-validation}

\item{is.binary.features}{the indicator whether the features are binary}

\item{k}{the number of features retained}

\item{subset.sel}{the indice of training cases used to select features}

\item{theta0}{the left bound of the uniform distribution as the prior of 
              "theta". And the right bound is set  to be "1-theta0" automatically}

\item{no.theta}{the parameter in Simpson rule used to evaluate  the integration
w.r.t. "theta".  The integrant is evaluated at 2*(no.theta)+1 points. }

\item{alpha.shape}{the shape parameter of the inverse Gamma, which is the prior
distribution of  "alpha"} 

\item{alpha.rate}{the rate parameter of the inverse Gamma, as above}

\item{no.alpha}{ the number of "alpha"'s used in approximating the integral with 
                 respect to "alpha", using mid-point rule}
		
\item{correct}{the indicator whether the correction method shall be applied}

\item{no.theta.adj}{the parameter of Simpson rule used to evaluate the 
                    integration with respect to "theta" in evaluating the 
		    adjustment factor
		    }
}

\value{

\item{prediction}{ a matrix showing the detailed prediction result: 
                   the 1st column being the true responses, 
		   the 2nd being the predicted responses, 
		   the 3rd being the predictive probabilities of class 1 
		   and the 4th being the indicator whether wrong prediction 
		   is made.}

\item{aml}{ the average minus log probabilities }

\item{error.rate}{ the ratio of wrong prediction }

\item{mse}{ the average square error of the predictive probabilities }

\item{summary.pred}{ tabular display of the predictive probabilities 
                     and the actual fraction of class 1. }

\item{alpha.prior.adj.post}{ a matrix showing the detailed information 
                             about the "alpha"'s, 
			     the 1st column being the values of "alpha"'s, 
			     the 2nd being the adjustment factor, i.e. 
			     probability that feature is discarded by 
			     the cutoff used in the feature selection, 
			     the 3rd being the log of the 2nd column times 
			     the numbers of discarded features, 
			     the 4th being the posterior probabilities}

\item{features.selected}{The features selected using correlation criterion}
}

\seealso{\code{\link{gendata.bayes}}}

\examples{

#generate a dataset 
d <- gendata.bayes(100,100,500,500,1000,400)

#do prediction with correction applied
pred.d.cor <- Predict.bayes(d$test,d$train,TRUE,10,,0,20,0.5,5,20,TRUE,40)

#do prediction without correction applied
pred.d.uncor <- Predict.bayes(d$test,d$train,TRUE,10,,0,20,0.5,5,20,FALSE,40)

#do 5-fold cross-validation on the training data with correction applied
cv.dtr.cor <- cv.bayes(d$train,TRUE,5,10,0,20,0.5,5,20,TRUE,40)

}

\references{ L.Li, J.Zhang, and R.M.Neal (2006) A Method for Avoiding Bias from
Feature Selection With Application to Naive Bayes Classification Models, Technical
Report 0705, Department of Statistics, University of Toronto}

\keyword{classif}
