\name{bootsPLS-package}
\alias{bootsPLS-package}
\docType{package}
\title{
bootsPLS: bootstrap subsamplings of sPLS-DA for classification and signature identification.}
\description{
This package helps for the discovery of new signatures by performing splsda on numerous random subsamplings; it provides tools to select only the most stable variables.
}

\details{
\tabular{ll}{
Package: \tab bootsPLS \cr
Type: \tab Package \cr
Version: \tab 1.0.0 \cr
Date: \tab 20-03-2015\cr
License: \tab GPL-3 \cr
}

The package implements the methodology described in Rohart \emph{et al.} (2015) for identifying genes that differentiate Mesemchymal Stromal Cells from other cell types. The method is applicable to any classification problem with more than 2 classes. It relies on bootstrap subsamplings of sPLS-DA and provides tools to select the most stable variables (defined as the ones consistently selected over the bootstrap subsamplings) and to predict the class of test samples.\cr\cr

Three major functions:\cr
*\code{\link{bootsPLS}} performs a \code{\link[mixOmics]{splsda}} of the \code{mixOmics} package on several random subsamplings. It records the selected variables on each replication. The \code{\link{compile.bootsPLS.object}} function is to be used if several calls to \code{\link{bootsPLS}} are made.  \cr
*\code{\link{fit.model}} fits a constraint spls (see \code{\link{spls.hybrid}}) on the most stable variables. The variables can be automatically selected by setting auto.tune=TRUE.\cr
*\code{\link{prediction}} predicts the class of test samples and compute Confidence Interval (CI) of the prediction values. \cr

}


\author{
Florian Rohart [aut,cre] \cr
Maintainer: florian.rohart@gmail.com
}


\references{
Rohart \emph{et al.} (2015). A Molecular Classification of Human Mesenchymal Stromal Cells\cr
Le Cao, K.-A., Boitard, S. and Besse, P. (2011). Sparse PLS Discriminant Analysis: biologically relevant feature selection and graphical displays for multiclass problems. \emph{BMC Bioinformatics} \bold{12}:253.

}


\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y


boot1=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)
boot2=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5,showProgress=FALSE)
boot3=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5,showProgress=FALSE)

# construct a list of bootsPLS object
bootsPLS.object=list(boot1,boot2,boot3)

# compile the outputs in one bootsPLS object
boot=compile.bootsPLS.object(bootsPLS.object)

# fit the model
fit=fit.model(boot,auto.tune=TRUE)
fit=fit.model(boot,ncomp=2)
plotIndiv(fit,ind.names=FALSE,col=as.numeric(Y))

# prediction and Confidence Interval
pred=prediction(fit,X.test=X,CI=TRUE)

head(pred$Y.hat.test[,,"comp.1"])
lapply(pred$out.CI$CI$'comp.1',head)

}
}
