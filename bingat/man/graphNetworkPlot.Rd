\name{graphNetworkPlot}
\alias{graphNetworkPlot}
\title{Graph Network Plots}

\description{This function plots the connections between nodes in a single subject.}

\usage{graphNetworkPlot(data, type, main = "Network Plot", labels, groupCounts, groupLabels)}

\arguments{
  \item{data}{A vector of a single graph.}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{main}{The title for the plot.}
  \item{labels}{A vector which contains the names for each node.}
  \item{groupCounts}{A vector which contains the number of nodes in each group of nodes.}
  \item{groupLabels}{A vector which contains the names for each group of nodes.}
}

\value{A plot displaying the connections between the nodes.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	main <- "Brain Connections"
	gc <- c(5, 5, 4, 6)
	gl <- c("Grp1", "Grp2", "Grp3", "Grp4")
	
	graphNetworkPlot(braingraphs[,1], "adjMatrix", main, groupCounts=gc, groupLabels=gl)
}
