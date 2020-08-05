\name{getNumNodes}
\alias{getNumNodes}
\title{Get the Number of Nodes in a Graph}

\description{This function will return the number of nodes for a given of graph.}

\usage{getNumNodes(data, type)}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
}

\value{The number of individual nodes in the given graph.}

\examples{
	data(braingraphs)
	
	brainnodes <- getNumNodes(braingraphs, "adjMatrix")
	brainnodes
}
