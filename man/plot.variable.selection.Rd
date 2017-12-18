\name{plot.variable.selection}
\alias{plot.variable.selection}
\docType{package}
\title{
Plot the results of the testing procedure to determine the number of variables to select}
\description{
Plot the results of the testing procedure to determine the number of variables to select on each component}
\usage{
\method{plot}{variable.selection}(x,pch,col,type,lty,...)
}

\arguments{
\item{x}{Object of class "variable.selection".}
\item{pch}{A vector of length ncomp (argument of \code{\link{bootsPLS}}). See \code{\link{par}} for details.}
\item{col}{A vector of length ncomp (argument of \code{\link{bootsPLS}}). See \code{\link{par}} for details.}
\item{type}{A vector of length ncomp (argument of \code{\link{bootsPLS}}). See \code{\link{par}} for details.}
\item{lty}{A vector of length ncomp (argument of \code{\link{bootsPLS}}). See \code{\link{par}} for details.}
\item{\dots}{not used}
}

\details{This function plots the p-value of the one sided t-test with respect to the frequency of each variable. See \code{\link{variable.selection}} and the reference below for more details on the multiple testing procedure.}

\references{
Rohart \emph{et al.} (2016). A Molecular Classification of Human Mesenchymal Stromal Cells. PeerJ, DOI 10.7717/peerj.1845\cr
}

\seealso{\code{\link{bootsPLS}}, \code{\link{variable.selection}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y


boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)
out=variable.selection(boot,ncomp=2,limit=c(40,40))
plot(out)


}
}
