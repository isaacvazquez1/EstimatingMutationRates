\name{ML_BivariateJackknife_CI}
\alias{ML_BivariateJackknife_CI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Jackknife (1-\eqn{\alpha})100\% confidence interval for maximum likelihood estimation of mutation rates with mutator strength fixed
%%  ~~function to do ... ~~
}
\description{The maximum likelihood estimation of mutation rates procedure is replicated \code{Replicates} times in parallel to obtain a jackknife (1-\eqn{\alpha})100\% confidence interval. Each replicate is performed with a random dilution factor whose lower and upper bounds are \code{dilution_lowerbound} and \code{dilution_upperbound}, respectively. The \code{ML_Jackknife_CI} function detects the number of available cores and \code{FreeCores} denotes the number of cores to keep unused for computer estability. Lastly, the mutation rate from mutator to mutant cells is regarded as the product of the mutation rate from wildtype to mutant cells and the scaling factor \code{strength}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ML_BivariateJackknife_CI(Replicates, alpha, FreeCores, File, population_size, resample_size, GC, dilution_lowerbound, dilution_upperbound, discretization, WT2R_bounds, WT2M_bounds, strength)
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
  \item{File}{A .txt file where the mutant cells and the mutant mutator cells from the growth cycle \eqn{j} appear on the \eqn{2j-1}-th and \eqn{2j}-th column, where \eqn{j=1, 2,\ldots, } \code{GC}.}
\item{population_size}{Final population size of the cell culture.}
%%     ~~Describe \code{x} here~~
\item{resampling}{\code{True/False} value to consider or not a resampling.}
\item{resample_size}{Number of cultures to work with, if \code{resampling=True}.}
\item{GC}{Growth cycle to work with.}
\item{dilution_lowerbound}{Dilution lower bound
%%     ~~Describe \code{dilution_lowerbound} here~~
}
  \item{dilution_upperbound}{Dilution upper bound
%%     ~~Describe \code{dilution_upperbound} here~~
}
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

FILENAME="QM-Mutants.txt"
ML_BivariateJackknife_CI(Replicates=500, alpha=0.05, FreeCores=5, File=FILENAME, population_size=2^33, resample_size=50, GC=1, dilution_lowerbound=0.00005, dilution_upperbound=2^20/2^33, discretization=25, WT2R_bounds=c(-8, -5), WT2M_bounds=c(-7, -3), strength=500)

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
