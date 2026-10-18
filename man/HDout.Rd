\name{HDout}
\alias{HDout}
\title{
Outlier detection in high-dimensional data
}
\description{
Distance-based outlier detection for high dimension, low sample size data (Ahn, Lee and Lee, 2018). The observation farthest from the rest is flagged as a potential outlier and removed, one at a time, up to \code{Mout} observations. Each step is displayed by plotting the observed distances against the distances expected without outliers, estimated by a parametric bootstrap. 
}
\usage{
HDout(data, method = 1:3, Mout = 4, Nsim = 30, plotY = TRUE)
}
\arguments{
  \item{data}{
\code{data} is a p by n data matrix (p > n).
}
  \item{method}{
choice of distance: \code{1} = centroid distance (\code{\link{dist_cen}}), \code{2} = ridge Mahalanobis (rM) distance (\code{\link{dist_mah}}), \code{3} = maximal data piling (MDP) distance (\code{\link{dist_mdp}}). The default \code{1:3} runs all three.
}
  \item{Mout}{
maximum number of outliers to test. The default is 4.
}
  \item{Nsim}{
number of bootstrap simulations. The default is 30; 10 is faster and 100 is more accurate.
}
  \item{plotY}{
if \code{TRUE} (default), the results are displayed in a plot.
} % item
} % argument
\details{
At each step (\code{Nout} = 1, ..., \code{Mout}), the leave-one-out distance of each remaining observation j is the distance between the single-member group \{j\} and the group of all other observations, computed with \code{\link{dist_cen}}, \code{\link{dist_mah}} or \code{\link{dist_mdp}}. For the rM distance, the ridge correction constant is \eqn{\sqrt{\log(p)/(m-1)}}, where m is the number of remaining observations. The observation with the largest distance is flagged as a potential outlier and removed.

The bootstrap samples are drawn from a normal distribution with the sample mean and covariance of the data without the flagged observation(s); when p >= n, a ridge term \eqn{\sqrt{\log(p)/m}\, I_p} is added to the covariance. For each of the \code{Nsim} samples, the sorted leave-one-out distances are computed and then averaged.

The plot has one row per method and one column per step. Each panel plots the observed sorted distances (y-axis) against the average simulated sorted distances (x-axis); the flagged observation is shown as a filled point labelled with its column index. A flagged point lying far above the pattern of the others suggests an outlier.
}
\value{
A list (returned invisibly) with the following elements:
\itemize{
\item \code{dj_c}: sorted centroid distances at the last iteration (\code{Mout}).
\item \code{dj_rM}: sorted rM distances at the last iteration.
\item \code{dj_MDP}: sorted MDP distances at the last iteration.
\item \code{outliers}: a list of the column indices flagged at steps 1, ..., \code{Mout}, one vector per method.
}
Distances are 0 for methods that were not run.
}
\references{
Ahn, J., Lee, M. H. and Lee, J. A. (2018). Distance-based outlier detection for high dimension, low sample size data, Journal of Applied Statistics. doi:10.1080/02664763.2018.1452901
}
\author{
Jung Ae Lee <jungae.lee@umassmed.edu>
}
\seealso{
\code{\link{dist_cen}}, \code{\link{dist_mah}}, \code{\link{dist_mdp}}
}
\examples{
# simulated data example: p = 1000 variables, n = 19 observations
data(simdata)
dim(simdata)   # 1000 by 19

y = simdata    # p by n matrix (p >= n)

# first 200 variables, small Mout and Nsim keep the example fast;
# HDout(y) with the defaults (Mout = 4, Nsim = 30) takes longer
out = HDout(y[1:200, ], Mout = 2, Nsim = 5)
out$outliers   # observations flagged at each step, per method

\donttest{
# leukemia data: first 200 genes, transposed to genes (p) by patients (n);
# centroid distance only, without the plot
data(leukemia)
HDout(t(leukemia$X[, 1:200]), method = 1, Mout = 2, Nsim = 5, plotY = FALSE)$outliers
}
}
