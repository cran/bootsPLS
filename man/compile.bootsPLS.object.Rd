\name{compile.bootsPLS.object}
\alias{compile.bootsPLS.object}
\docType{package}
\title{
Combine several bootsPLS objects into one}
\description{
Combine several bootsPLS objects into one}
\usage{
compile.bootsPLS.object(bootsPLS.list,path,pattern,file,save.file)

}

\arguments{
\item{bootsPLS.list}{A list of bootsPLS object, as obtained from multiple calls to \code{\link{bootsPLS}}}
\item{path}{Only used if \code{bootsPLS.list} is missing. A path to the Rdata files containing the bootsPLS outputs to combine. Will be passed in the \code{list.files} function.}
\item{pattern}{Only used if \code{bootsPLS.list} is missing. A commun pattern to the Rdata files to combine. Will be passed in the \code{list.files} function.}
\item{file}{Only used if \code{bootsPLS.list} is missing. Vector of Rdata files to combine. Will be loaded as \code{path/file[i]}}
\item{save.file}{Optional. Full path of the Rdata file to be saved with the combined objects.}
}


\details{This function works with either a list of \code{bootsPLS} object, a \code{path} and a \code{pattern} or a \code{path} and a \code{file}. The outputs can be saved into a Rdata file. See Examples below.}

\value{A 'bootsPLS' object is returned for which \code{\link{plot}}, \code{\link{fit.model}} and \code{\link{prediction}} are available.
    See \code{\link{bootsPLS}} for the outputs.

}



\seealso{\code{\link{plot.bootsPLS}}, \code{\link{fit.model}}, \code{\link{prediction}}, \code{\link{bootsPLS}}}

\examples{
\dontrun{
data(MSC)
X=MSC$X
Y=MSC$Y
dim(X)
table(Y)

# perform several bootsPLS analysis
boot1=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5)
boot2=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5,showProgress=FALSE)
boot3=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5,showProgress=FALSE)

# construct a list of bootsPLS object
bootsPLS.object=list(boot1,boot2,boot3)

# compile the outputs in one bootsPLS object
boot=compile.bootsPLS.object(bootsPLS.object)


# =======================================
# all can work from the saved file:
# =======================================
#convenient if used on a cluster

# perform several bootsPLS analysis and save the outputs
boot1=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5,save.file="MSC1.Rdata")
boot2=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5,showProgress=FALSE,save.file="MSC2.Rdata")
boot3=bootsPLS(X=X,Y=Y,ncomp=3,many=5,kCV=5,showProgress=FALSE,save.file="MSC3.Rdata")

# compile the outputs in one bootsPLS object
boot=compile.bootsPLS.object(path=paste(getwd(),"/",sep=""),
    pattern="MSC",save.file="MSC.all.Rdata")

# or
boot=compile.bootsPLS.object(path=paste(getwd(),"/",sep=""),
    file=c("MSC1.Rdata","MSC2.Rdata","MSC3.Rdata"),save.file="MSC.all.Rdata")


}
}
