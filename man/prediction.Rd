\name{prediction}
\alias{prediction}
\docType{package}
\title{
prediction}
\description{
prediction}
\usage{
prediction (object,X,Y,signature,ncomp,X.test,CI,many,
                subsampling.matrix,ratio,level.CI,save.file)
}

\arguments{

\item{object}{a `spls.constraint' object, as one resulting from \code{\link{fit.model}}. If object is missing: X, Y, signature are needed.}

\item{X}{Only used if \code{object} is missing. Input matrix of dimension n * p; each row is an observation vector.}
\item{Y}{Only used if \code{object} is missing. Factor with at least q>2 levels; Response variable of length n, indicating the sample types.}
\item{signature}{Only used if \code{object} is missing. A list containing which variables are to be kept on each component.}
\item{ncomp}{Only used if \code{object} is missing. How many component do you want to include in the sPLS-DA analysis?}
\item{X.test}{Test matrix.}
\item{CI}{logical. If TRUE, the confidence interval are calculated.}
\item{many}{How many subsamplings do you want to do? Default is 100}
\item{subsampling.matrix}{Optional matrix of \code{many} columns. Gives the samples to subsample as an internal learning set.}
\item{ratio}{Number between 0 and 1. It is the proportion of the n samples that are put aside and considered as an internal testing set. The (1-ratio)*n samples are used as a training set.}
\item{level.CI}{A 1- \code{level.CI}\% confidence interval is calculated.}
\item{save.file}{Save the outputs of the functions in \code{save.file.Rdata}.}

}


\details{
This function can work with a \code{spls.constraint} object or with the input data (X, Y, signature). See examples below to see the difference in use.

}

\value{


\item{CI}{A (1- \code{level.CI})\% confidence interval is returned for each samples in \code{X.test}}
\item{Y.hat.test}{A four dimensional array. The two first dimensions are an estimation of the dummy matrix obtained from Y (size n * number of sample types). The third dimension is relative to the number of components \code{ncomp}. The fourth dimension concerns the number of subsamplings.}
\item{ClassifResult}{A 5-dimensional array. The two first dimensions consists in the confusion matrix. The third dimension is relative to the number of components \code{ncomp}. The fourth dimension concerns the number of subsamplings. The fifth and last dimension is relative to the different distances "max.dist", "centroids.dist" and "mahalanobis.dist". }
\item{loadings.X}{A 3-dimensional array. Loadings vector of X, for each component and each subsampling.}
\item{prediction.X}{A 4-dimensional array of size n*many*ncomp*3. Gives the prediction for the chosen  \code{method} of all the samples, either in the internal learning set or the internal testing set. The last dimension is relative to the different distances "max.dist", "centroids.dist" and "mahalanobis.dist".}
\item{prediction.X.test}{A 4-dimensional array of size nrow(X.test)*many*ncomp*3. Gives the prediction for the chosen  \code{method} of all the test samples in X.test. The last dimension is relative to the different distances "max.dist", "centroids.dist" and "mahalanobis.dist".}
\item{learning.sample}{Matrix of size n*many. Gives the samples that have been used in the internal training set over the \code{many} replications. These samples have the value 1, the others 0. }
\item{coeff}{A list of means.X, sigma.X, means.Y and sigma.Y. Means and variances for the variables of X and the columns of the dummy matrix obtained from Y, each row is a subsampling.}
\item{data}{A list of the input data X, Y and of ind.kept.X, which is a list containing the variables kept on each component.}


}

\seealso{\code{\link{fit.model}}, \code{\link{CI.prediction}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y


# with a bootsPLS object
boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)
fit=fit.model(boot,ncomp=3)

# with a spls.constraint object and without CI
pred=prediction(fit,X.test=X)

lapply(pred$predicted.test,head)

# with a spls.constraint object and with CI
pred.CI=prediction(fit,X.test=X,CI=TRUE)
lapply(pred.CI$out.CI$CI$'comp.1',head)
lapply(pred.CI$out.CI$CI$'comp.2',head)
lapply(pred.CI$out.CI$CI$'comp.3',head)

# without a spls.constraint object. X,Y and signature are needed
# the results should be similar
#(not the same because of the random subsamplings,
# exactly the same if subsampling.matrix is an input)
signature=fit$data$signature
pred=prediction(X=X,Y=Y,signature=signature,X.test=X)

pred2=prediction(X=X,Y=Y,signature=signature,X.test=X,CI=TRUE)
lapply(pred2$out.CI$CI$'comp.1',head)
lapply(pred2$out.CI$CI$'comp.2',head)
lapply(pred2$out.CI$CI$'comp.3',head)

}
}
