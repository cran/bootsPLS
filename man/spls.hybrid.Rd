\name{spls.hybrid}
\alias{spls.hybrid}
\docType{package}
\title{
spls.hybrid, midway between PLS and sPLS}
\description{
Performs a constraint sPLS on the first PLS-components and a sPLS on the last components}
\usage{

spls.hybrid(X,
            Y,
            ncomp = 2,
            mode = c("regression", "canonical", "invariant", "classic"),
            max.iter = 500,
            tol = 1e-06,
            keepX.constraint,
            keepY.constraint,
            keepX,
            keepY,
            near.zero.var = FALSE)
}

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{numeric vector or matrix of responses (for multi-response models).
\code{NA}s are allowed.}
\item{ncomp}{the number of components to include in the model (see Details).
Default is 2.}
\item{mode}{character string. What type of algorithm to use.
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}. See Details.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{tol}{a positive real, the tolerance used in the iterative algorithm.}
\item{keepX.constraint}{A list containing which variables of X are to be kept on each of the first PLS-components.}
\item{keepY.constraint}{A list containing which variables of Y are to be kept on each of the first PLS-components.}
\item{keepX}{number of \eqn{X} variables kept in the model on the last components.}
\item{keepY}{number of \eqn{Y} variables kept in the model on the last components.}
\item{near.zero.var}{boolean, see the internal \code{\link[mixOmics]{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations.}

}


\details{The \code{spls.hybrid} function allows you to compute a constraint spls on the first components and a spls on the last components. Note that the only condition on \code{keepX.constraint} and \code{keepX} is that the sum of both length is ncomp; likewise for the ones relative to Y.}

\value{A 'spls.hybrid' object is returned. The object is a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{Y}{the centered and standardized original response vector or matrix.}
\item{ncomp}{the number of components included in the model.}
\item{mode}{the algorithm used to fit the model.}
\item{keepX.constraint}{A list of length \code{ncomp} containing which variables of X are to be kept on each component.}
\item{keepY.constraint}{A list of length \code{ncomp} containing which variables of Y are to be kept on each component.}
\item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
\item{variates}{list containing the variates.}
\item{loadings}{list containing the estimated loadings for the \eqn{X} and
\eqn{Y} variates.}
\item{names}{list containing the names to be used for individuals and variables.}
\item{nzv}{list containing the zero- or near-zero predictors information, for \code{X} and \code{Y}.}
\item{coeff}{A list of means.X, sigma.X, means.Y and sigma.Y. Means and variances for the variables of \code{X} and the columns of \code{Y.mat}.}
}


\references{
Rohart \emph{et al.} (2016). A Molecular Classification of Human Mesenchymal Stromal Cells. PeerJ, DOI 10.7717/peerj.1845\cr
}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y

Y.mat=unmap(Y)

n=nrow(X)
p=ncol(X)


keepX.constraint=list(sample(1:p,7),sample(1:p,15))
fit=spls.hybrid(X,Y=Y.mat,ncomp=4,keepX.constraint=keepX.constraint,keepX=c(5,10))

plotIndiv(fit,ind.names=FALSE)
}
}
