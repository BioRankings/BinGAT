\name{genAlgConsensus}
\alias{genAlgConsensus}
\alias{gaConsensus}
\title{Find Edges Separating Two Groups using Multiple Genetic Algorithm's (GA) Consensus}
\description{
	GA-Mantel is a fully multivariate method that uses a genetic algorithm to search over 
	possible edge subsets using the Mantel correlation as the scoring measure for assessing 
	the quality of any given edge subset.
}

\usage{
	genAlgConsensus(data, covars, consensus = .5, numRuns = 10, 
		parallel = FALSE, cores = 3, ...)
}

\arguments{
  \item{data}{A matrix of edges(rows) for each sample(columns).}
  \item{covars}{A  matrix of covariates(columns) for each sample(rows).}
  \item{consensus}{The required fraction (0, 1] of solutions containing an edge in order to keep it.}
  \item{numRuns}{
  		Number of runs to do.  In practice the number of runs needed varies based on data set size
  		and the GA parameters set.
  }
  \item{parallel}{When this is 'TRUE' it allows for parallel calculation of the bootstraps. Requires the package \code{doParallel}.}
  \item{cores}{The number of parallel processes to run if parallel is 'TRUE'.}
  \item{...}{Other arguments for the GA function see \link[bingat]{genAlg}}
}

\value{
	A list containing

	\item{solutions}{The best solution from each run.}
	\item{consSol}{The consensus solution.}
	\item{selectedIndex}{The selected edges by row number.}
}

\details{
	Use a GA consensus approach to find edges that separate subjects based on group membership 
	or set of covariates if you cannot run the GA long enough to get a final solution.
}

\author{Sharina Carter, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
	\dontrun{
		data(braingraphs)
		
		### Set covars to just be group membership
		covars <- matrix(c(rep(0, 19), rep(1, 19)))
		
		### We use low numbers for speed. The exact numbers to use depend
		### on the data being used, but generally the higher iters and popSize 
		### the longer it will take to run.  earlyStop is then used to stop the
		### run early if the results aren't improving.
		iters <- 500
		popSize <- 200
		earlyStop <- 250
		numRuns <- 3
		
		gaRes <- genAlgConsensus(braingraphs, covars, .5, numRuns, FALSE, 3, 
				iters, popSize, earlyStop)
	}
}