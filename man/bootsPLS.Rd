\name{bootsPLS}
\alias{bootsPLS}
\docType{package}
\title{
Performs replications of the splsda on random subsamplings of the data}
\description{
Performs replications of the splsda on random subsamplings of the data}
\usage{bootsPLS(X,Y,near.zero.var,many,ncomp,
            method = c("max.dist", "centroids.dist", "mahalanobis.dist"),
            save.file,ratio,kCV,grid,showProgress=TRUE)
}

\arguments{
\item{X}{Input matrix of dimension n * p; each row is an observation vector.}
\item{Y}{Factor with at least q>2 levels.}
\item{near.zero.var}{Logical. If TRUE, a pre-screening step is performed to remove predictors with near-zero variance. See \code{\link[mixOmics]{nearZeroVar}}.}
\item{many}{How many replications of the sPLS-DA analysis are to be done?}
\item{ncomp}{How many component are to be included in the sPLS-DA analysis?}
\item{method}{Indicates the distance that is used to classify the samples. One of "max.dist", "centroids.dist", "mahalanobis.dist". Default is "max.dist"}
\item{save.file}{If the outputs are to be saved, this argument allows you to do it at the end of each replication. A full path is expected. Convenient if you run this function on a cluster and it is killed before completion, e.g. due to a too short requested time.}
\item{ratio}{Number between 0 and 1. It is the proportion of the n samples that are put aside and considered as an internal testing set. The (1-ratio)*n samples are used as a training set and the \code{kCV} fold cross validation is performed on them. Default is 0.3}
\item{kCV}{Number of fold for the cross validation. Default is 10.}
\item{grid}{List of length ncomp. grid[[i]] contains a grid of value for the tuning of the \code{keepX} parameter of the sPLS-DA on component i. See \code{\link[mixOmics]{spls}} for more details on \code{keepX}. Default is \code{grid[[i]]=1:min(40,ncol(X))}.  }
\item{showProgress}{Logical. If TRUE, shows the progress of the algorithm: which one of the \code{many} replications and which of the \code{kCV} fold of the cross validation are currently done. It also gives a list of which variables are selected on each component.}
}


\details{Performs replication of the \code{\link[mixOmics]{splsda}} on random subsamplings of the data and record which variables are selected on which subsamplings. It also gives a confusion matrix for each component and for each subsamplings.}

\value{A 'bootsPLS' object is returned for which \code{plot}, \code{fit.model} and \code{prediction} are available.


\item{ClassifResult}{A 4-dimensional array. The two first dimensions consists in the confusion matrix. The third dimension is relative to the number of components \code{ncomp}. The fourth dimension concerns the number of replication \code{many}.  }
\item{loadings.X}{A 3-dimensional array. Loadings vector of X, for each component and each replication.}
\item{selection.variable}{A 3-dimensional array. Gives the selected variables for each component and each replication. It is obtained by replacing each non zero value in \code{loadings.X} by 1.}
\item{frequency}{A matrix of size ncomp*p. Gives the frequency of selection for each variable on each component. It is obtained as a mean over the third dimension of \code{selection.variable}}
\item{nbr.var}{Matrix of size many*ncomp. Gives the number of variables that have been selected on each component for each replication.}
\item{learning.sample}{Matrix of size n*many. Gives the samples that have been used in the internal training set over the \code{many} replications. These samples have the value 1, the others 0. }
\item{prediction}{A 3-dimensional array of size n*many*ncomp. Gives the prediction for the chosen  \code{method} of all the samples, either in the learning set or the testing set.}
\item{data}{A list of the input data X, Y and of the method used to classify the sample ("max.dist", "centroids.dist" or "mahalanobis.dist").}

}

\references{
Rohart \emph{et al.} (2015). A Molecular Classification of Human Mesenchymal Stromal Cells\cr
}

\seealso{\code{\link[mixOmics]{splsda}}, \code{\link{plot.bootsPLS}}, \code{\link{fit.model}}, \code{\link{prediction}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y
dim(X)
table(Y)


boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)


# saving the outputs in a Rdata file
save.file=paste(getwd(),"/MSC.",Sys.getpid(),".Rdata",sep="")
boot=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5,save.file=save.file)

}
}
