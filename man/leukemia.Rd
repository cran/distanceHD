\name{leukemia}
\alias{leukemia}
\docType{data}
\title{
Gene expression data from Golub et al. (1999)
}
\description{
Gene expression data (3051 genes and 38 tumor mRNA samples) from the leukemia microarray study of Golub et al. (1999). )
}
\usage{data(leukemia)}
\format{
A list with the following elements:
\itemize{
\item \code{X}: a (38 by 3051) matrix giving the expression levels of 3051 genes for 38 leukemia patients. Each row corresponds to a patient, each column to a gene.
\item \code{Y}: a numeric vector of length 38 giving the cancer class of each patient.
\item \code{gene.names}: a matrix containing the names of the 3051 genes for the gene expression matrix X. The three columns correspond to the gene index, ID, and Name, respectively.
}
}
\source{
The data are described in Golub et al. (1999) and can be freely downloaded from http://www.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?paper_id=43. 
}
\references{
S. Dudoit, J. Fridlyand and T. P. Speed (2002). Comparison of discrimination methods for the classification of tumors using gene expression data, Journal of the American Statistical Association 97, 77-87.

Golub et al. (1999). Molecular classification of cancer: class discovery and class prediction by gene expression monitoring, Science 286, 531-537.
}
\seealso{plsgenomics::leukemia}
\examples{
data(leukemia)

# how many samples and genes?
dim(leukemia$X)

# how many samples of class 1 and 2, respectively?
table(leukemia$Y)
}
\keyword{datasets}
