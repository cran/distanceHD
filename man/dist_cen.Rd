\name{dist_cen}
\alias{dist_cen}
\alias{dist_centroid}
\title{
Centroid distance between two groups
}
\description{
Calculate the Centroid (Euclidean) distance between two groups in a high-dimensional space (d > n), with support for a single-member cluster. It also works in low-dimensional settings.
}
\usage{
dist_cen(x, group)
}
\arguments{
  \item{x}{
\code{x} is n by d matrix
}
  \item{group}{
\code{group} is a binary group label with the length of n1 and n2
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
x = leukemia$X          # 38 by 3051 genes

# apply the function
dist_cen(x, group)  # 25.4
}
