\name{calcDistance}
\alias{calcDistance}
\title{Calculate the Distance Between Vectors}

\description{This function calculates the distance between two vectors.}

\usage{calcDistance(x, y, type = "", method = "hamming")}

\arguments{
  \item{x, y}{Vectors of the same length that contain 1's and 0's.}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt). See 'Details'}
  \item{method}{The distance metric to use, currently only "hamming" is supported.}
}

\details{
	If the \code{type = "adjMatrix"} is used, the value will be divided by 2 to account for duplicate comparisons.
	Otherwise the \code{type} does not affect the output.
}

\value{A single number indicating the distance between the two input vectors.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	dist <- calcDistance(braingraphs[,1], braingraphs[,2], "adjMatrix")
	dist
}
