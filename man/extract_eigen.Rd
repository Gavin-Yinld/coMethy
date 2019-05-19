\name{extract_eigen}
\alias{extract_eigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
extract eigen-loci from co-methylation module
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
extract eigen-loci from co-methylation module

}
\usage{
extract_eigen(csm.ml, all_label, number_of_eig)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csm.ml}{
%%     ~~Describe \code{csm.ml} here~~
The methylation profile with each row denotes one loci and each column denotes a sample
}
  \item{all_label}{
%%     ~~Describe \code{all_label} here~~
A vector contains the module labels for all the loci
}
  \item{number_of_eig}{
%%     ~~Describe \code{number_of_eig} here~~
Total number of eigen-loci need to be extracted
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
#eigen_loci <- extract_eigen(methy_data=module$profile,all_label=module$module_id,number_of_eig=100,plot=T)
}

