\name{evaluate_by_loss}

\alias{evaluate_by_loss}

\title{calculating the total loss of prediction results}

\description{Calculates the average loss of the predictive probabilities by assuming the ratio of the loss of assigning 0 to 1 to the loss of assigning 1 to 0 is \code{ratio.loss}. Note that the prediction threshold is no longer 0.5. "1" will be assigned if the predictive probability is greater than \code{1/(1+ratio.loss)}.
}


\usage{
evaluate_by_loss ( y.true, pred.prob, ratio.loss=10)
}

\arguments{

\item{y.true}{a vector containing the true response.}

\item{pred.prob}{a vector containing the predictive probabilities.}

\item{ratio.loss}{ratio of loss1to0 to loss0to1.}

}

\value{

\item{loss}{the average loss.}

\item{sd}{the standard deviation of the average loss.}

}

\keyword{classif}
