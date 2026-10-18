\name{simdata}
\alias{simdata}
\docType{data}
\title{
Simulated high-dimensional data for outlier detection
}
\description{
A simulated high dimension, low sample size data set with 1000 variables and 19 observations, used to illustrate the outlier detection procedure \code{\link{HDout}}. 
}
\usage{data(simdata)}
\format{
A numeric matrix of dimension 1000 by 19 (p by n, p >= n). Each row corresponds to a variable and each column to an observation, which is the input format of \code{\link{HDout}}.
}
\seealso{
\code{\link{HDout}}
}
\examples{
data(simdata)

# how many variables (p) and observations (n)?
dim(simdata)    # 1000 x 19

\donttest{
y = simdata     # p by n matrix (p >= n)

out = HDout(y[1:200, ], Mout = 2, Nsim = 5)   # first 200 variables, for speed
out$outliers   # observations flagged at each step, per method
}
}
\keyword{datasets}
