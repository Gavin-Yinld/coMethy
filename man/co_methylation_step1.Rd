\name{co_methylation_step1}
\alias{co_methylation_step1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The first step of co-methylation analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The first step of the co-methylation analysis, including kmeans analysis to group pCSM loci into three clusters, i.e. hypo/mid/hyper-methylation cluster, and, for each kmeans cluster, the network topology analysis function in WGCNA package is called to pick the soft-thresholding power.
}
\usage{
co_methylation_step1(csm_ml_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csm_ml_matrix}{
%%     ~~Describe \code{x} here~~
methylation profile of pCSM loci in each sample
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list the following two components:
\item{profile }{the methylation profile of pCSM loci}
\item{modult_id }{the label tells that which co-methylation module the pCSM loci belong to}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library("coMethy")
file=paste(system.file(package="coMethy"),"extdata/co_methy.test.data.txt",sep='/')
meth_data <- read.table(file,sep='\t',header=T,row.names=1)
kmeans_cluster <- co_methylation_step1(meth_data)
}

