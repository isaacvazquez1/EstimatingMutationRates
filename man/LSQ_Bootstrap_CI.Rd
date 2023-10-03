\name{LSQ_Bootstrap_CI}
\alias{LSQ_Bootstrap_CI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap (1-\eqn{\alpha})100\% confidence interval for direct estimation of mutation rates
%%  ~~function to do ... ~~
}
\description{The direct estimation of mutation rates procedure is replicated \code{Replicates} times in parallel to obtain a bootstrap (1-\eqn{\alpha})100\% confidence interval. The \code{LSQ_Bootstrap_CI} function detects the number of available cores and \code{FreeCores} denotes the number of cores to keep unused for computer estability.
}
\usage{
LSQ_Bootstrap_CI(Replicates, alpha, FreeCores, File, population_size, GC, dilution, discretization, WT2R_bounds, WT2M_bounds, M2R_bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Replicates}{Number of replicates to obtaing the bootstrap (1-\eqn{\alpha})100\% confidence interval.
%%     ~~Describe \code{Replicates} here~~
}
  \item{alpha}{Value between 0 and 1.
%%     ~~Describe \code{alpha} here~~
}
  \item{FreeCores}{Cores to keep unused.
%%     ~~Describe \code{FreeCores} here~~
}
  \item{File}{A .txt file where the mutant cells and the mutant mutator cells from the growth cycle \eqn{j} appear on the \eqn{2j-1}-th and \eqn{2j}-th column, where \eqn{j=1, 2,\ldots, } \code{GC}.
%%     ~~Describe \code{File} here~~
}

  \item{population_size}{Final population size of the cell culture.}
%%     ~~Describe \code{x} here~~

\item{GC}{Growth cycle to work with.}
\item{dilution}{Dilution factor.}
\item{discretization}{Number of pieces, minus one, each mutation rate sample space will be splitted into.}
\item{WT2R_bounds}{A two-entry numerical vector.}
\item{WT2M_bounds}{A two-entry numerical vector.}
\item{M2R_bounds}{A two-entry numerical vector.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Bootstrap (1-\eqn{\alpha})100\% confidence interval for mutation rate from wildtype to mutant cells, from wildtype to mutator cells, and from mutator to mutant cells.
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


FILENAME="GM-Mutants.txt"
LSQ_Bootstrap_CI(Replicates=500, alpha=0.05, FreeCores=5, File=FILENAME, population_size=2^33, GC=1, dilution=0.0005, discretization=25, WT2R_bounds=c(-9.5, -5.5), WT2M_bounds=c(-8, -4), M2R_bounds=c(-7, -3))


## The function is currently defined as
function (Replicates, alpha, FreeCores, File, population_size,
    GC, dilution, discretization, WT2R_bounds, WT2M_bounds, M2R_bounds)
{
    library(parallel)
    cl <- makeCluster(detectCores() - FreeCores)
    cl
    library(doParallel)
    registerDoParallel(cl)
    rep = Replicates
    output <- foreach(k = 1:rep) \%dopar\% {
        EstimatingMutationRates::LSQ_Estimation(File, population_size,
            GC, dilution, discretization, WT2R_bounds, WT2M_bounds,
            M2R_bounds)
    }
    result_1 = c()
    for (j in 1:rep) {
        result_1 = rbind(result_1, output[[j]])
    }
    estimate_WT2R = mean(result_1[, 1])
    estimate_WT2M = mean(result_1[, 2])
    estimate_M2R = mean(result_1[, 3])
    confidence_intervals = matrix(0, 3, 3)
    confidence_intervals[1, 1] = estimate_WT2R
    confidence_intervals[2, 1] = estimate_WT2M
    confidence_intervals[3, 1] = estimate_M2R
    alpha
    for (i in 1:ncol(result_1)) {
        aux = result_1[, i]
        aux = sort(aux)
        lower = round(alpha/2 * (nrow(result_1) + 1))
        upper = round((1 - alpha/2) * (nrow(result_1) + 1))
        if (lower == 0) {
            lower = 1
        }
        if (upper > nrow(result_1)) {
            upper = nrow(result_1)
        }
        confidence_intervals[i, 2] = aux[lower]
        confidence_intervals[i, 3] = aux[upper]
    }
    colnames(confidence_intervals) = c("Estimate", "Lower bound",
        "Upper bound")
    rownames(confidence_intervals) = c(paste("Wildtype to Mutant | GC",
        GC), paste("Wildtype to Mutator | CG", GC), paste("Mutator to Mutant mutator  | GC",
        GC))
    return(confidence_intervals)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
