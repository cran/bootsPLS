\name{MSC}
\alias{MSC}
\docType{data}
\title{
Mesenchymal Stem Cells data}
\description{
The data is a subset of the data contained in Rohart et al. (2015). It contains 150 samples, either Mesenchymal Stem Cells (MSC) or Non-MSC, and 200 genes.}

\usage{
data(MSC)
}


\format{
A list containing the following components:
\describe{
\item{\code{X}}{data matrix with 150 rows and 200 columns. Each row represents
an experimental sample, and each column a single gene.}
\item{\code{Y}}{a vector containing the sample type (MSC or Non-MSC).}
}
}

\details{ This data is a subset of a bigger study; the samples and genes have been randomly selected for this package. All samples have been background corrected, log2 transformed and YuGene normalised.
}

\references{
Rohart \emph{et al.} (2015). A Molecular Classification of Human Mesenchymal Stromal Cells\cr
Kim-Anh Le Cao, Florian Rohart, Leo McHugh, Othmar Korn, Christine A. Wells. YuGene: A simple approach to scale gene expression data derived from different platforms for integrated analyses. Genomics. http://dx.doi.org/10.1016/j.ygeno.2014.03.001.\cr
}
