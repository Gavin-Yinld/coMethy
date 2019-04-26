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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (csm.ml, all_label, number_of_eig) 
{
    eig_num_per_module <- function(label, total_number = 10000) {
        freq <- table(label)
        base_line <- floor(total_number/length(freq))
        min_cluster_num <- min(freq)
        pick_number <- 0
        result <- data.frame(module_id = names(freq), number = 0)
        rownames(result) <- names(freq)
        picked_num <- sum(result$number)
        while (picked_num < total_number) {
            if (base_line > min_cluster_num) {
                index <- which(freq < base_line)
                result[names(freq)[index], 2] <- freq[index]
                picked_num <- sum(result$number)
                freq <- freq[-index]
                base_line <- floor((total_number - picked_num)/length(freq))
                min_cluster_num <- min(freq)
            }
            else {
                result[names(freq), 2] <- base_line
                picked_num <- sum(result$number)
                if (picked_num < total_number) {
                  rest <- total_number - picked_num
                  index <- order(freq, decreasing = TRUE)[1:rest]
                  result[names(freq)[index], 2] <- result[names(freq)[index], 
                    2] + 1
                  picked_num <- sum(result$number)
                }
            }
        }
        if (picked_num > total_number) {
            rest <- picked_num - total_number
            index <- order(freq, decreasing = TRUE)[1:rest]
            result[index, 2] <- result[index, 2] - 1
            picked_num <- sum(result$number)
        }
        return(result)
    }
    number_per_module <- eig_num_per_module(label = all_label, 
        total_number = number_of_eig)
    nmf.input = NULL
    nmf.input.label <- NULL
    for (label in unique(all_label)) {
        temp.data <- all.data[which(all_label == label), ]
        p <- prcomp(t(temp.data), scale = FALSE)
        temp.res <- temp.data[order(abs(p$rotation[, 1]), decreasing = T)[1:number_per_module[label, 
            2]], ]
        nmf.input <- rbind(nmf.input, temp.res)
        nmf.input.label <- c(nmf.input.label, rep(label, nrow(temp.res)))
    }
    return(list(methy_prof = nmf.input, nmf.input.label = nmf.input.label))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
