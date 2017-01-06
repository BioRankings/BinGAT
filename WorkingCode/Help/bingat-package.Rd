\name{bingat-package}
\alias{bingat-package}
\alias{bingat}
\docType{package}
\title{Binary Graph Analysis Tools}

\description{
	Tools for analyzing binary graphs, including calculating the MLE of a set of binary graphs, 
	comparing MLE of sets of graphs, regression analysis on sets of graphs, using genetic algorithm 
	to identify nodes and edges separating sets of graphs, and generating random binary graphs 
	sampled from the Gibbs distribution.
}

\details{
	The following are the types of binary graphs that are accepted:
	\enumerate{
		\item adjMatrix: An entire binary adjacency matrix as a single vector
		\item adjMatrixLT: The upper or lower triangle of a binary adjacency matrix as a single vector
		\item diag: The diagonal vector on a binary adjacency matrix
	}	
}

\references{
	\enumerate{
		\item Stat Med. 2015 Nov 25. doi: 10.1002/sim.6757.
				Gibbs distribution for statistical analysis of graphical data with a sample application to fcMRI brain images.
				La Rosa PS1,2, Brooks TL1, Deych E1, Shands B1,3, Prior F4, Larson-Prior LJ4,5, Shannon WD1,3.
	}
}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\keyword{ package }




