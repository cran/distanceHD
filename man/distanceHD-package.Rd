\name{distanceHD-package}
\alias{distanceHD-package}
\alias{distanceHD}
\docType{package}
\title{
Distance Metrics for High-Dimensional Clustering
}
\description{
We provide three distance metrics for measuring the separation between two clusters in high-dimensional spaces. The first metric is the centroid distance, which calculates the Euclidean distance between the centers of the two groups. The second is a ridge Mahalanobis distance, which incorporates a ridge correction constant, alpha, to ensure that the covariance matrix is invertible. The third metric is the maximal data piling distance, which computes the orthogonal distance between the affine spaces spanned by each class. These three distances are asymptotically interconnected and are applicable in tasks such as discrimination, clustering, and outlier detection in high-dimensional settings.
}

\author{
Jung Ae Lee <jungae.lee@umassmed.edu>; Jeongyoun Ahn <jyahn@kaist.ac.kr>

Maintainer: Jung Ae Lee <jungaeleeb@gmail.com>
}
\references{
1. Ahn J, Marron JS (2010). The maximal data piling direction for discrimination. Biometrika, 97(1):254-259.

2. Ahn J, Lee MH, Yoon YJ (2012). Clustering high dimension, low sample size data using the maximal data piling distance. Statistica Sinica, 22(2):443-464.

3. Ahn J, Lee MH, Lee JA (2019). Distance-based outlier detection for high dimension, low sample size data. Journal of Applied Statistics.46(1):13-29.
}
\keyword{ centroid distance }
\keyword{ ridge Mahalanobis distance }
\keyword{ maximal data piling distance }
