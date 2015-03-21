\name{plot.bootsPLS}
\alias{plot.bootsPLS}
\docType{package}
\title{
Plot the frequency of selection of all variables for all the PLS-component.}
\description{
Plot the frequency of selection of all variables for all the PLS-component.}
\usage{
\method{plot}{bootsPLS}(x,light,pch,col,legend.position,title,...)
}

\arguments{
\item{x}{Object of class `bootsPLS'.}
\item{light}{Between 0 and 1. A lot of variables usually have a low frequency, this argument allows you to discard those from the plot. Only the variables that have a frequency higher than \code{light} on at least one component are plotted.}
\item{pch}{A vector of length ncomp (argument of \code{\link{bootsPLS}}). See \code{\link{par}} for details.}
\item{col}{A vector of length ncomp (argument of \code{\link{bootsPLS}}). See \code{\link{par}} for details.}
\item{legend.position}{Location of the legend, see \code{\link{legend}} for more details. If set to FALSE, no legend is plotted.}
\item{title}{title of the plot}
\item{\dots}{not used}
}

\details{Plot the frequency of selection of all variables depending on the component. By default, a legend and a title are added to the plot, those can be changed/removed by changing the respective inputs.}


\seealso{\code{\link{bootsPLS}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y
dim(X)
table(Y)


boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)

plot(boot)
plot(boot,light=0.3,legend.position=FALSE)
}
}
