\name{fit.model}
\alias{fit.model}
\docType{package}
\title{
Create a spls.constraint object by fitting a constraint spls on a bootsPLS object}
\description{
Create a spls.constraint object by fitting a constraint spls with the \code{spls.hybrid} function on a `bootsPLS' object}
\usage{
fit.model(object,auto.tune,X,Y,ncomp,signature,alpha,limit,showProgress=TRUE)

}

\arguments{
\item{object}{a `bootsPLS' object', as obtained from \code{\link{bootsPLS}}.}
\item{auto.tune}{Logical. If TRUE, tune the optimal number of component (ncomp) and which variables to select on each component (signature). It only works with \code{object}}
\item{X}{Input matrix of dimension n * p; each row is an observation vector.}
\item{Y}{Factor with at least q>2 levels.}

\item{ncomp}{How many component are to be included in the sPLS-DA analysis?}
\item{signature}{A list containing which variables to keep on each component.}
\item{alpha}{Level of the test.}
\item{limit}{Vector of maximal number of genes to include on each component.}
\item{showProgress}{Logical. If TRUE, shows the progress of the algorithm.}
}



\details{This function fit a \code{spls.hybrid} on the variables included in \code{signature}, which can be an input or internally calculated by setting \code{auto.tune=TRUE}. If \code{object} is given as an input, (X, Y) are ignored. If \code{auto.tune=TRUE}, \code{ncomp, signature} are ignored.}

\value{A 'spls.constraint' object is returned for which \code{\link[mixOmics]{plotIndiv}} is available.

The outputs are the ones from \code{\link{spls.hybrid}}, plus
\item{data}{A list of the input data X, Y, Y.mat (dummy matrix) and of \code{signature}, which is a list containing the variables kept on each component. If a tuning occurs, either for the number of components or for the variables or both, outputs $component.selection and $variable.selection are available; see \code{\link{component.selection}} and \code{\link{variable.selection}}.}


}

\seealso{\code{\link{prediction}}, \code{\link{CI.prediction}}, \code{\link[mixOmics]{plotIndiv}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y


boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)

# with a bootsPLS object and auto.tune=TRUE
fit=fit.model(boot,auto.tune=TRUE)
plot(fit$component.selection)
plot(fit$variable.selection)

# with a bootsPLS object and ncomp=2
fit=fit.model(boot,ncomp=2)

# with a bootsPLS object and ncomp/signature as input
signature=fit$data$signature
fit=fit.model(boot,ncomp=2,signature=signature)

# with no bootsPLS object
fit=fit.model(X=X,Y=Y,ncomp=2,signature=signature)# bootsPLS object

plotIndiv(fit,ind.names=FALSE, legend=TRUE)

}
}
