\name{QuantileM_GrowthCycles}
\alias{QuantileM_GrowthCycles}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quantile function model for \eqn{n} growth cycles}
\description{The quantile function method allows us to simulate phenotypic dynamics under the mutator phenotype hypothesis.

Let \code{mu_WR}, \code{mu_WM}, and \code{mu_MR} be the mutation rate from wiltype cells to mutant cells, the mutation rate from wildtype cells to mutator cells, and the mutation rate from mutator cells to mutant mutator cells, respectively. Suppose a cell culture is started with a quantity \code{initial_pop} of wildtype cells.

After the cell culture has reached the population size of \code{final_pop} cells, i.e., after completing the first growth cycle, a random sample of size \code{bottle_neck} cells will be taken and used to start a new cell culture which has to reach a final population size of \code{final_pop} cells. This process will be repeated as many times as \code{n_cycles} indicates.

Outputs:
\itemize{
\item{Mutations leading to mutant cells (from wildtype cells) at growth cycle \eqn{j}}
\item{Mutant cells at growth cycle \eqn{j}}
\item{Mutations leading to mutator cells (from wildtype cells) at growth cycle \eqn{j}}
\item{Mutator cells at growth cycle \eqn{j}}
\item{Mutations leading to mutant mutator cells (from wildtype cells) at growth cycle \eqn{j}}
\item{Mutant mutator cells at growth cycle \eqn{j}}
}
for each \eqn{j=1, \ldots, } \code{n_cycles}

Remark.

\code{initial_pop} \eqn{\in \mathbb{N}.}

\code{final_pop} \eqn{\in \mathbb{N}.}

\code{bottle_neck} \eqn{\in \mathbb{N}.}

\code{initial_pop} \eqn{<} \code{bottle_neck} \eqn{<} \code{final_pop}.

\eqn{0 \leq } \code{mu_WR} \eqn{ \leq 1}.

\eqn{0 \leq } \code{mu_WM} \eqn{ \leq 1}.

\eqn{0 \leq } \code{mu_MR} \eqn{ \leq 1}.
}
\usage{
QuantileM_GrowthCycles(initial_pop, final_pop, mu_WR, mu_WM, mu_MR, n_cycles, bottle_neck)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{initial_pop}{Initial inoculum of the cell culture.}
  \item{final_pop}{Final population size of the cell culture.}
  \item{bottle_neck}{Number of cells in the random sample taken from the cell culture to start the next gowth cycle.}
  \item{n_cycles}{Number of growth cycles}
  \item{mu_WR}{Mutation rate from wildtype cells to mutant cells.}
  \item{mu_WM}{Mutation rate from wildtype cells to mutator cells.}
  \item{mu_MR}{Mutation rate from mutator cells to mutant mutator cells.}
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
QuantileM_GrowthCycles(initial_pop=1, final_pop=10^9, mu_WR=10^-7, mu_WM=10^-5, mu_MR=10^-5, n_cycles=3, bottle_neck=10^{3})



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
