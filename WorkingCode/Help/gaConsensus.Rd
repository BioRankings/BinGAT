\name{gaConsensus}
\alias{gaConsensus}
\title{Genetic Algorithm Consensus Solution}

\description{Runs a Genetic Algorithm to find nodes or edges that separate groups}

\usage{
	gaConsensus(data, groups, iters = 10, nRuns = 1, popSize = 200, 
		method = "manhattan", parallel = FALSE, cores = 3)
}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{groups}{A vector of 0/1s used to denote the group membership of columns in the data.}
  \item{iters}{The number of iterations (ideally, 1000 or more).}
  \item{nRuns}{The number of times to repeat the GA to get a consensus.}
  \item{popSize}{The population size.}
  \item{method}{The distance metric to use.}
  \item{parallel}{TRUE or FALSE depending on whether the analysis will be parallelized for speed.}
  \item{cores}{The number of cores to use for parallelization. Ignored if parallel = FALSE.}
}

\seealso{
  \code{\link[genalg]{rbga.bin}}.
}

\value{A list containing all the solutions from the last run(s) of the GA and their correlation values}

\author{Berkley Shands, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	grps <- c(rep(0, 19), rep(1, 19))
	iters <- 1 ### This is set low for speed
	nRuns <- 1 ### This is set low for speed
	
	consensus <- gaConsensus(braingraphs, grps, iters, nRuns) 
	consensus$corr[1:5]
	consensus$solutions[1:2,]
}
