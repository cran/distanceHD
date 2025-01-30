\name{dist_mah}
\alias{dist_mah}
\alias{dist_ma}
\title{
ridge Mahalanobis distance between two groups
}
\description{
Calculate the Mahalanobis distance between two groups in a high-dimensional space (d > n), using a ridge correction on the covariance matrix to ensure invertibility. The method also supports a single-member cluster and works in low-dimensional settings without a ridge correction.
}
\usage{
dist_mah(x, group, alpha)
}
\arguments{
  \item{x}{
\code{x} is n by d matrix
}
  \item{group}{
\code{group} is a binary group label with the length of n1 and n2
}
\item{alpha}{
\code{alpha} is a positive numeric value representing the ridge correction constant. If not specified, the default value is set to \eqn{\sqrt{\log(d)/n}}.
} % item
} % argument
\value{
A numeric value of distance
}
\author{
Jung Ae Lee <jungae.lee@umassmed.edu>
}

\examples{
data(leukemia)
group = leukemia$Y      # 38 patients status with a value of 1 or 2
x = leukemia$X          # 38 x 3051 genes

# apply the function
dist_mah(x, group)   # 26.8 

# default alpha
d = 3051; n = 38
alpha = sqrt(log(d)/n)
dist_mah(x, group, alpha)   # 26.8
}

