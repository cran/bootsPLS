\name{plot.predictCI}
\alias{plot.predictCI}
\docType{package}
\title{
Plot confidence Intervals}
\description{
Display the confidence intervals for the prediction values of a specific level of the outcome}
\usage{
\method{plot}{predictCI}(x,
ncomp=1, level=1, las=2, col, title, name.var = TRUE, abline=TRUE,...)
}

\arguments{

\item{x}{a `predict.CI' object, as one resulting from \code{\link{CI.prediction}}, or from \code{\link{prediction}$out.CI}.}
\item{ncomp}{Component to be plotted}
\item{level}{One of the levels of the outcome to plot}
\item{las}{las argument for the x-axis labels}
\item{col}{Color for the confidence Intervals. Can be a single value or a different value for each sample.}
\item{title}{title of the plot}
\item{name.var}{If TRUE, add the name of the variables on the x-axis. A vector of length the number of variables can be supplied for customized names.}
\item{abline}{logical. Add two horizontal lines: 0.5 and 1/the number of levels of the outcome. Note that if the outcome only has two levels, only one line is plotted.}

\item{\dots}{not used}

}


\details{Plot the Confidence Intervals for a specific component and a specific level. In the special case of the "max.dist" being used (in the \code{bootsPLS} function), a line at 0.5 and 1/number of levels is added to the plot with \code{abline=TRUE}. These lines highlight the predicted class of the samples. For example if the level 1 out of 2 levels is displayed: any sample above 0.5 is predicted as level 1, any samples below 0.5 is predicted as level 2, and samples with a confidence interval overlaying 0.5 are not a clear cut. Another example if the level 1 out of 3 levels is displayed: any sample above 0.5 is predicted as level 1, any samples below 1/3 is predicted as not level 1, and samples in between 1/3 and 0.5 can't be decided from this graph only (you can display other levels, or find the predicted class in your object).
    See example below.

    In the case of another distance being used, the lines should not be added to the graph \code{abline=FALSE}.

    }

\seealso{\code{\link{CI.prediction}}, \code{\link{prediction}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y


# with a bootsPLS object
boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)
fit=fit.model(boot,ncomp=3)

CI=CI.prediction(fit,X.test=X)
plot(CI)

pred=prediction(fit,X.test=X, CI=TRUE)
plot(pred$out.CI, ncomp=2)
# we color each confidence interval by the predicted class
plot(pred$out.CI, ncomp=2,
    col = color.mixo(factor(pred$predicted.test$"max.dist"[,2]))) #second component


# because bootsPLS was used with dist="max.dist" (by default)
# and because there are only two levels in the outcome,
# everything above the 0.5 line is predicted as MSC

par(mfrow=c(3,1))
plot(pred$out.CI, ncomp=1, level="MSC")
plot(pred$out.CI, ncomp=2, level="MSC")
plot(pred$out.CI, ncomp=3, level="MSC")


#we can do the same things for the second level (Non-MSC)
par(mfrow=c(3,1))
plot(pred$out.CI, ncomp=1, level="Non-MSC")
plot(pred$out.CI, ncomp=2, level="Non-MSC")
plot(pred$out.CI, ncomp=3, level="Non-MSC")

}
}
