\name{ML_Bivariate}
\alias{ML_Bivariate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum likelihood estimation for mutator strength fixed
%%  ~~function to do ... ~~
}
\description{Given a growth cycle \code{GC} where \eqn{c} cell cultures where plated and each culture contains \code{population_size} cells, a dilution process is performed according a binomial sampling, this is, a random sample of size one is drawn from a Binomial distribution whose parameters are \eqn{n}, the number on mutant cells in each culture in the growth cycle, and \eqn{p=}\code{dilution}, the dilution factor.

From the above framework a resample of size \code{resample_size} with replacement is taken to consider a smaller quanty of experiments. The maximum likelihood estimation methods search over a parameter space defined by the user. Vectors \code{WT2R_bounds}, \code{WT2M_bounds} will be containing the lower and upper bound on log\eqn{_{10}} scale of the mutation rates on their first and second entry, respectively. Each vector is representing an interval which be splitted into \code{discretization}\eqn{+1} pieces to be a computationally tractable problem. Lastly, the mutation rate from mutator to mutant cells is regarded as the product of the mutation rate from wildtype to mutant cells and the scaling factor \code{strength}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ML_Bivariate(File, population_size, resampling, resample_size, GC, dilution, discretization, WT2R_bounds, WT2M_bounds, strength)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  \item{strength}{A positive factor.
%%     ~~Describe \code{strength} here~~
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

FILENAME="GM-Mutants.txt"
ML_Bivariate(File=FILENAME, population_size=2^33, resampling = F, resample_size=0, GC=1, dilution=0.0001, discretization=25, WT2R_bounds=c(-8, -5), WT2M_bounds=c(-7, -3), strength= 1000)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (File, population_size, resampling, resample_size, GC,
    dilution, discretization, WT2R_bounds, WT2M_bounds, strength)
{
    binomial_sampling_GC_j <- function(td, j, ps1) {
        a = c()
        for (i in 1:nrow(td)) {
            a = c(a, rbinom(n = 1, size = (td[, 2 * j - 1] +
                td[, 2 * j])[i], prob = ps1))
        }
        return(a)
    }
    pmf_coeff <- function(k, WT2R, WT2M, M2R, pop, N, ps1, strength) {
        sum = 0
        for (n in 0:(N - 1)) {
            aux_char = exp(((0+2i) * pi * n)/N)
            aux = (-1 * (0+2i) * pi * n * k)/N
            sum = sum + pgs(aux_char, p = ps1, Am = 10^{
                WT2R
            } * pop, m1 = 10^{
                WT2M
            }, Amm = strength * 10^{
                WT2R
            } * pop) * exp(aux)
        }
        return(Re(sum/N))
    }
    discretization_fun <- function(range, n) {
        mu_min = range[1]
        mu_max = range[2]
        a = seq(0, n)
        result = mu_min * (mu_max/mu_min)^{
            a/n
        }
        return(result)
    }
    ll_fun <- function(data, WT2R, WT2M, M2R, pop_size, N, ps1,
        strength) {
        sum(log(pmf_coeff(data + 1, WT2R, WT2M, M2R, pop_size,
            N, ps1, strength)))
    }
    td <- read.table(File, header = F, sep = "")
    data <- binomial_sampling_GC_j(td, j = GC, p = dilution)
    if (resampling == T) {
        data = sample(data, resample_size, replace = T)
    }
    N = max(data)
    N
    paramters_WT2R = WT2R_bounds
    WT2R = discretization_fun(paramters_WT2R, discretization)
    paramters_WT2M = WT2M_bounds
    WT2M = discretization_fun(paramters_WT2M, discretization)
    max = ll_fun(data, WT2R = WT2R[1], WT2M = WT2M[1], M2R = WT2R[1],
        population_size, N, dilution, strength)
    powers = c()
    output = c()
    for (i in WT2R) {
        for (j in WT2M) {
            aux = ll_fun(data, WT2R = i, WT2M = j, M2R = i, population_size,
                N, dilution, strength)
            output = rbind(output, c(10^{
                i
            }, 10^{
                j
            }, strength * 10^{
                i
            }, aux))
            if (max < aux) {
                max = aux
                powers = c(10^{
                  i
                }, 10^{
                  j
                }, strength * 10^{
                  i
                })
            }
        }
    }
    max
    names(powers) = c("Wildtype to Mutant", "Wildtype to Mutator",
        "Mutator to Mutant mutator")
    df = as.data.frame(output)
    df_confidence <- df[df[, 4] >= max - qchisq(0.95, 2)/2, ]
    CI = matrix(0, 3, 2)
    row.names(CI) = c("Wildtype to Mutant", "Wildtype to Mutator",
        "Mutator to Mutant mutator")
    colnames(CI) = c("Lower bound", "Upper bound")
    CI[1, 1] = min(df_confidence$V1)
    CI[1, 2] = max(df_confidence$V1)
    CI[2, 1] = min(df_confidence$V2)
    CI[2, 2] = max(df_confidence$V2)
    CI[3, 1] = min(df_confidence$V3)
    CI[3, 2] = max(df_confidence$V3)
    out = list(powers, CI, max, output)
    return(out)
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
