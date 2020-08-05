\name{estMLE}
\alias{estMLE}
\title{Estimate the MLE Parameters}

\description{This function estimates the MLE parameters g-star and tau for a given set of graphs.}

\usage{estMLE(data, type)}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
}

\details{
	Essentially this function calls both \code{estGStar} and \code{estTau} and returns the results. \cr
}

\value{A list containing g-star and tau named gstar and tau respectively.}

\examples{
	data(braingraphs)

	brainmle <- estMLE(braingraphs, "adjMatrix") 
	brainmle
}
