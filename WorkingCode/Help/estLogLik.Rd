\name{estLogLik}
\alias{estLogLik}
\title{Estimate the Log Likelihood Value}

\description{This function estimates log likelihood value for a given graph.}

\usage{estLogLik(data, type, gstar, tau, g = NULL)}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{gstar}{A single vector to estimate the likelihood for.}
  \item{tau}{A single value used in estimating the likelihood.}
  \item{g}{Deprecated. Replaced with gstar for clarity.}
}

\value{The log-likelihood value of the data.}

\examples{
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	brainll <- estLogLik(braingraphs, "adjMatrix", braingstar, braintau)
	brainll
}
