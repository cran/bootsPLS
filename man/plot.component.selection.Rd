\name{plot.component.selection}
\alias{plot.component.selection}
\docType{package}
\title{
Plot the results of the testing procedure to determine the number of component to select}
\description{
Plot the results of the testing procedure to determine the number of component to select}
\usage{
\method{plot}{component.selection}(x,pch,col,type,lty,...)
}

\arguments{
\item{x}{Object of class "component.selection".}
\item{pch}{ See \code{\link{par}} for details.}
\item{col}{ See \code{\link{par}} for details.}
\item{type}{See \code{\link{par}} for details.}
\item{lty}{See \code{\link{par}} for details.}
\item{\dots}{not used}
}

\details{This function plots the p-value of the one sided t-test with respect to the number of components. See \code{\link{component.selection}} and the reference below for more details on the multiple testing procedure}

\references{
Rohart \emph{et al.} (2015). A Molecular Classification of Human Mesenchymal Stromal Cells\cr
}

\seealso{\code{\link{bootsPLS}}, \code{\link{component.selection}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y


boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)
comp=component.selection(boot)

plot(comp)

}
}
