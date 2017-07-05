\name{genAlg}
\alias{genAlg}
\title{Find Edges Separating Two Groups using Genetic Algorithm (GA)}
\description{
	GA-Mantel is a fully multivariate method that uses a genetic algorithm to search over 
	possible edge subsets using the Mantel correlation as the scoring measure for assessing 
	the quality of any given edge subset.
}

\usage{
	genAlg(data, covars, iters = 50, popSize = 200, earlyStop = 0, 
		dataDist = "manhattan", covarDist = "gower", verbose = FALSE, 
		plot = TRUE, minSolLen = NULL, maxSolLen = NULL)
}

\arguments{
  \item{data}{A matrix of edges(rows) for each sample(columns).}
  \item{covars}{A  matrix of covariates(columns) for each sample(rows).}
  \item{iters}{The number of times to run through the GA.}
  \item{popSize}{The number of solutions to test on each iteration.}
  \item{earlyStop}{
  		The number of consecutive iterations without finding a better solution before stopping 
  		regardless of the number of iterations remaining. A value of '0' will prevent early 
  		stopping.
  		}
  \item{dataDist}{The distance metric to use for the data. This can only be "manhattan" for now.}
  \item{covarDist}{The distance metric to use for the covariates. Either "euclidean" or "gower".}
  \item{verbose}{While 'TRUE' the current status of the GA will be printed periodically.}
  \item{plot}{A boolean to plot the progress of the scoring statistics by iteration.}
  \item{minSolLen}{The minimum number of columns to select.}
  \item{maxSolLen}{The maximum number of columns to select.}

}

\value{
	A list containing

	\item{scoreSumm}{
		A matrix summarizing the score of the population.  This can be used to figure out if the ga has 
		come to a final solution or not. This data is also plotted if plot is 'TRUE'.
	}
	\item{solutions}{The final set of solutions, sorted with the highest scoring first.}
	\item{scores}{The scores for the final set of solutions.}
	\item{time}{How long in seconds the ga took to run.}
	\item{selected}{The selected edges by name.}
	\item{nonSelected}{The edges that were NOT selected by name.}
	\item{selectedIndex}{The selected edges by row number.}
}

\details{
	Use a GA approach to find edges that separate subjects based on group membership or set of covariates.
	
	The data and covariates should be normalized BEFORE use with this function because of distance
	functions.
	
	This function uses modified code from the rbga function in the genalg package. \link[genalg]{rbga}
	
	Because the GA looks at combinations and uses the raw data, edges with a small difference 
	may be selected and large differences may not be.
	
	The distance calculations use the vegdist package. \link[vegan]{vegdist}
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
		
		gaRes <- genAlg(braingraphs, covars, iters, popSize, earlyStop)
	}
}