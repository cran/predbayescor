\name{evaluate_by_loss}

\alias{evaluate_by_loss}

\title{calculating the total loss of prediction results}

\description{Calculates the average loss of predictions based on threshold with \code{threshold}. Note that this threshold has 1-1 mapping with the ratio of the loss of assigning 0 to 1 to the loss of assigning 1 to 0: threshold=1-1/(1+ratio). 
}


\usage{
evaluate_by_loss ( y.true, pred.prob, threshold=0.5)
}

\arguments{

\item{y.true}{a vector containing the true response.}

\item{pred.prob}{a vector containing the predictive probabilities.}

\item{threshold}{When predictive probability is greater than \code{threshold}, the response is predicted as 1.}

}

\value{

\item{loss}{the average loss, with attrib "sd" storing the estimate of the standard error of this loss.}

}

\keyword{classif}
