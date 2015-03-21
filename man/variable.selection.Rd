\name{variable.selection}
\alias{variable.selection}
\docType{package}
\title{
Tune the number of variables on each component}
\description{
Performs a testing procedure to choose the number of variable on each component}
\usage{
variable.selection(object,ncomp,alpha,limit)

}


\arguments{
\item{object}{a `bootsPLS' object', as obtained from \code{\link{bootsPLS}}}
\item{ncomp}{How many component?}
\item{alpha}{Level of the test}
\item{limit}{Vector of maximal number of genes to include on each component}

}


\details{The testing procedure ranks the variables by decreasing frequency in \code{object$frequency}, for each component. Random subsamplings are constructed, a \code{spls.hybrid} is fitted on the internal learning set and a prediction is made on the internal test set.
The testing procedure evaluates the gain in classification accuracy over the random subsamplings when a new variable is added from a decreasing frequency. This is done by on-sided t-test of level \code{alpha}. See the reference below for more details on the multiple testing procedure.

}

\value{A 'variable.selection' object is returned for which \code{plot} is available.

\item{pval}{pvalue obtained from the testing procedure.}
\item{opt}{Number of components chosen by the procedure}
\item{keepX.constraint}{Variables chosen on each of the \code{ncomp} components}
\item{object}{input `bootsPLS' object}
\item{alpha}{input level of the test}


}

\references{
Rohart \emph{et al.} (2015). A Molecular Classification of Human Mesenchymal Stromal Cells\cr
}

\seealso{\code{\link{spls.hybrid}}, \code{\link{plot}}, \code{\link{fit.model}}, \code{\link{prediction}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y
dim(X)
table(Y)

boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)
out=variable.selection(boot,ncomp=2)

out2=variable.selection(boot,ncomp=2,limit=c(40,40))

}
}
