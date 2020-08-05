\name{estTau}
\alias{estTau}
\title{Estimate Tau}

\description{This function estimates tau for a given set of graphs.}

\usage{estTau(data, type, gstar)}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{gstar}{A single columned data frame to be used as the g-star of the data set.}
}

\value{The tau value for the data based on g star.}

\examples{
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	braintau
}
