\name{getGibbsMixture}
\alias{getGibbsMixture}
\title{Group Splitter}

\description{This function splits the data into groups based on the Gibbs criteria.}

\usage{getGibbsMixture(data, type, desiredGroups, maxIter = 50, digits = 3)}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{desiredGroups}{The number of groups to test for.}
  \item{maxIter}{The maximum number of iterations to run searching for an optimal split.}
  \item{digits}{The number of digits to round internal values to when checking the stop criteria.}
}

\details{
	Generally this function is not used by itself but in conjunction with getLoglikeMixture.
}

\value{
	A list that contains information about the group splits.  The list contains the final weights,
	gstars and taus for every group, a boolean indicating convergence, the number of iterations
	it took, and the group for each graph.
}

\examples{
	data(braingraphs)

	braingm <- getGibbsMixture(braingraphs, "adjMatrix", 5) 
}
