\name{getNumEdges}
\alias{getNumEdges}
\title{Get the Number of Edges in a Graph}

\description{This function will return the number of edges for a given of graph.}

\usage{getNumEdges(nodes, type)}

\arguments{
  \item{nodes}{The number of individual nodes in a given graph.}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
}

\value{The number of edges between individual nodes in the given graph.}

\examples{	
	data(braingraphs)
	
	brainnodes <- getNumNodes(braingraphs, "adjMatrix")
	brainedges <- getNumEdges(brainnodes, "adjMatrix")
	brainedges
}
