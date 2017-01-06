\name{estGStar}
\alias{estGStar}
\title{Estimate G-Star}

\description{This function estimates the g-star graph for a given set of graphs.}

\usage{estGStar(data)}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
}

\value{A single vector that is the gstar is returned.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)

	braingstar <- estGStar(braingraphs) 
	braingstar[1:25]
}