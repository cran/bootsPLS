\name{component.selection}
\alias{component.selection}
\docType{package}
\title{
Tune the number of components}
\description{
Performs a multiple hypotheses testing procedure to choose the number of components in the splsda analysis}
\usage{
component.selection(object,alpha,showProgress=TRUE)

}

\arguments{
\item{object}{a `bootsPLS' object, as obtained from \code{\link{bootsPLS}}}
\item{alpha}{Level of the test}
\item{showProgress}{Logical. If TRUE, shows the progress of the algorithm.}

}


\details{The testing procedure evaluates the gain in classification accuracy when a new PLS-component is added. This is done by on-sided t-test of level \code{alpha} applied on the classification accuracy obtained in the `bootsPLS' object. See the reference below for more details on the multiple testing procedure.
}

\value{A 'component.selection' object is returned for which \code{plot} is available.

\item{pval}{pvalue obtained from the testing procedure.}
\item{opt}{Number of components chosen by the procedure}
\item{object}{input `bootsPLS' object}
\item{alpha}{input level of the test}


}

\references{
Rohart \emph{et al.} (2016). A Molecular Classification of Human Mesenchymal Stromal Cells. PeerJ, DOI 10.7717/peerj.1845\cr
}

\seealso{\code{\link{plot.component.selection}}, \code{\link{fit.model}}, \code{\link{prediction}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y
dim(X)
table(Y)

boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)

comp=component.selection(boot)
}
}
