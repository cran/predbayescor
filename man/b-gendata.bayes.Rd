\name{gendata.bayes}

\alias{gendata.bayes}

\title{Generate binary data with Bayesian naive Bayes Models}

\description{"gendata.bayes" generates data (both training and test data) with
Bayesian naive Bayes model. The prior distribution of "theta" is uniform(0,1). The value of "alpha" is given by argument alpha, which controls the the overall relationship between the response and the predictor variables.}

\usage{gendata.bayes(n0,n1,m0,m1,p,alpha)}

\arguments{
\item{n0}{the number of class 0 in training data}
\item{n1}{the number of class 1 in training data}
\item{m0}{the number of class 0 in test data}
\item{m1}{the number of class 1 in test data}
\item{p}{the number of features}
\item{alpha}{a parameter controlling the dependency between the features and the
             response}
}
\value{
\item{train}{the training data, with the row standing for the cases and the
first column being the response}
\item{test}{the test data, of the same format as "train"}
}

\seealso{
\code{\link{predict_bayes}}
}
\keyword{datagen}
