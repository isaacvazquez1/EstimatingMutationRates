\name{LSQ_Estimation}
\alias{LSQ_Estimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum likelihood estimation}
\description{Given a growth cycle \code{GC} where \eqn{c} cell cultures where plated and each culture contains \code{population_size} cells, a dilution process is performed according a binomial sampling, this is, a random sample of size one is drawn from a Binomial distribution whose parameters are \eqn{n}, the number on mutant cells in each culture in the growth cycle, and \eqn{p=}\code{dilution}, the dilution factor.

From the above framework a resample of size \code{resample_size} with replacement is taken to consider a smaller quanty of experiments. The maximum likelihood estimation methods search over a parameter space defined by the user. Vectors \code{WT2R_bounds}, \code{WT2M_bounds}, \code{M2R_bounds} will be containing the lower and upper bound on log\eqn{_{10}} scale of the mutation rates on their first and second entry, respectively. Each vector is representing an interval which be splitted into \code{discretization}\eqn{+1} pieces to be a computationally tractable problem.

Remark.
\code{population_size} \eqn{\in \mathbb{N}.}

\code{resample_size} \eqn{\in \mathbb{N}.}

\code{GC} \eqn{\in \mathbb{N}.}

\eqn{0\leq }\code{dilution} \eqn{\leq 1}

\code{discretization} \eqn{\in \mathbb{N}.}

\code{WT2R_bounds} \eqn{\in \mathbb{R}^2.}

\code{WT2M_bounds} \eqn{\in \mathbb{R}^2.}

\code{M2R_bounds} \eqn{\in \mathbb{R}^2.}

\code{population_size} \eqn{\geq} \code{resample_size}}
\usage{
ML_Estimation(File, population_size, resample_size, GC, dilution, discretization,
                     WT2R_bounds, WT2M_bounds, M2R_bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{File}{A .txt file where the mutant cells and the mutant mutator cells from the growth cycle \eqn{j} appear on the \eqn{2j-1}-th and \eqn{2j}-th column, where \eqn{j=1, 2,\ldots, } \code{GC}.}
\item{population_size}{Final population size of the cell culture.}
%%     ~~Describe \code{x} here~~
\item{resample_size}{Number of cultures to work with.}
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

#Generating 100 simultations
experiments=10
result=c()
for(i in 1:experiments){
  result=rbind(result, QuantileM_GrowthCycles(initial_pop=1, final_pop=2^33, mu_WR=10^-7,
  mu_WM=10^-5, mu_MR=10^-4, n_cycles=3, bottle_neck=2^20))
  }

#Saving the dataset
filename <- "QM - InitialPop 1, FinalPop 2^33, Bottleneck 2^20, Growth cycles 3, WT2R 10^-7, WT2M 10^-5, M2R 10^-4.txt"
write.table(result, filename, sep = " ",append = T,col.names = TRUE, row.names = FALSE)

#Reading the dataset
aux <- read.table(filename,header = T,sep=" ")

#Saving dataset of mutants of each growth cycle
Mutants=cbind(aux$Mutants_GC_1, aux$Mutant_mutators_GC_1, aux$Mutants_GC_2, aux$Mutant_mutators_GC_2,
      aux$Mutants_GC_3, aux$Mutant_mutators_GC_3)
FILENAME="QM-Mutants.txt"
write.table(Mutants, FILENAME, sep = " ",append = T,col.names = FALSE,row.names = FALSE)

#Maximum likelihood estimation
estimates=ML_Estimation(File=FILENAME, population_size=2^33, resample_size=50, GC=2, dilution=0.0001,
               discretization=25, WT2R_bounds=c(-8, -5), WT2M_bounds=c(-7, -3), M2R_bounds=c(-7, -3))

#Estimation of the mutation rates
mutation_rates=10^estimates
mutation_rates

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
