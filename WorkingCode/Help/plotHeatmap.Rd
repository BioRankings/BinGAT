\name{plotHeatmap}
\alias{plotHeatmap}
\title{Plot Heat Map}

\description{This function plots the connections between nodes in a single subject as a heat map.}

\usage{plotHeatmap(data, type, names, ...)}

\arguments{
  \item{data}{A vector of a single graph.}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{names}{A vector of names for labeling the nodes on the plot.}
  \item{...}{Arguments to be passed to the plot method.}
}

\value{A plot displaying the connections between the nodes as a heat map.}

\examples{
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	plotHeatmap(braingstar, "adjMatrix")
}
