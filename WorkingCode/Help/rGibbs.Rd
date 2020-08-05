\name{rGibbs}
\alias{rGibbs}
\title{Generate Random Data}

\description{Generate random data sampled from the Gibbs distribution.}

\usage{rGibbs(gstar, tau, type, numGraphs = 1)}

\arguments{
  \item{gstar}{G star vector.}
  \item{tau}{A single value that affects the dispersion of the generated data.}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{numGraphs}{The number of graphs to generate.}
}

\value{A data frame containing all the graphs generated.}

\examples{
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs)
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	randombraingraphs <- rGibbs(braingstar, braintau, "adjMatrix", 3) 
	randombraingraphs[1:5,]
}
