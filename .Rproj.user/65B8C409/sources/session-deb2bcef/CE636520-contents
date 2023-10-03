\name{ML_Jackknife_CI}
\alias{ML_Jackknife_CI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Jackknife (1-\eqn{\alpha})100\% confidence interval for maximum likelihood estimation of mutation rates
%%  ~~function to do ... ~~
}
\description{The maximum likelihood estimation of mutation rates procedure is replicated \code{Replicates} times in parallel to obtain a jackknife (1-\eqn{\alpha})100\% confidence interval. The \code{ML_Jackknife_CI} function detects the number of available cores and \code{FreeCores} denotes the number of cores to keep unused for computer estability.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ML_Jackknife_CI(Replicates, alpha, FreeCores, File, population_size, resample_size, GC, dilution_lowerbound, dilution_upperbound, discretization, WT2R_bounds, WT2M_bounds, M2R_bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Replicates}{Number of replicates to obtaing the jackknife (1-\eqn{\alpha})100\% confidence interval.
%%     ~~Describe \code{Replicates} here~~
}
  \item{alpha}{Value between 0 and 1.
%%     ~~Describe \code{alpha} here~~
}
  \item{FreeCores}{Cores to keep unused.
%%     ~~Describe \code{FreeCores} here~~
}
  \item{File}{A .txt file where the mutant cells and the mutant mutator cells from the growth cycle \eqn{j} appear on the \eqn{2j-1}-th and \eqn{2j}-th column, where \eqn{j=1, 2,\ldots, } \code{GC}.}
\item{population_size}{Final population size of the cell culture.}
%%     ~~Describe \code{x} here~~
\item{resampling}{\code{True/False} value to consider or not a resampling.}
\item{resample_size}{Number of cultures to work with, if \code{resampling=True}.}
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

FILENAME="QM-Mutants.txt"
ML_Jackknife_CI(Replicates=5, alpha=0.05, FreeCores=5, File=FILENAME, population_size=2^33, resample_size=50, GC=1, dilution_lowerbound=0.00005, dilution_upperbound=2^20/2^33, discretization=25, WT2R_bounds=c(-8, -5), WT2M_bounds=c(-7, -3), M2R_bounds=c(-6, -2))

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Replicates, alpha, FreeCores, File, population_size,
    resample_size, GC, dilution_lowerbound, dilution_upperbound,
    discretization, WT2R_bounds, WT2M_bounds, M2R_bounds)
{
    library(parallel)
    cl <- makeCluster(detectCores() - FreeCores)
    cl
    library(doParallel)
    registerDoParallel(cl)
    rep = Replicates
    output <- foreach(k = 1:rep) \%dopar\% {
        EstimatingMutationRates::ML_Estimation(File, population_size,
            resampling = T, resample_size, GC, dilution = runif(1,
                dilution_lowerbound, dilution_upperbound), discretization,
            WT2R_bounds, WT2M_bounds, M2R_bounds)
    }
    result_1 = c()
    for (j in 1:rep) {
        result_1 = rbind(result_1, output[[j]][[1]])
    }
    estimate_WT2R = mean(result_1[, 1])
    estimate_WT2M = mean(result_1[, 2])
    estimate_M2R = mean(result_1[, 3])
    confidence_intervals = matrix(0, 3, 3)
    confidence_intervals[1, 1] = estimate_WT2R
    confidence_intervals[2, 1] = estimate_WT2M
    confidence_intervals[3, 1] = estimate_M2R
    aux_WT2R = qt(p = alpha/2, df = rep - 1, lower.tail = F) *
        var(result_1[, 1])/sqrt(rep)
    aux_WT2M = qt(p = alpha/2, df = rep - 1, lower.tail = F) *
        var(result_1[, 2])/sqrt(rep)
    aux_M2R = qt(p = alpha/2, df = rep - 1, lower.tail = F) *
        var(result_1[, 3])/sqrt(rep)
    confidence_intervals[1, 2] = estimate_WT2R - aux_WT2R
    confidence_intervals[1, 3] = estimate_WT2R + aux_WT2R
    confidence_intervals[2, 2] = estimate_WT2M - aux_WT2M
    confidence_intervals[2, 3] = estimate_WT2M + aux_WT2M
    confidence_intervals[3, 2] = estimate_M2R - aux_M2R
    confidence_intervals[3, 3] = estimate_M2R + aux_M2R
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
