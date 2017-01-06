\name{braingraphs}
\alias{braingraphs}
\docType{data}
\title{Brain Graph Data Set}
\description{A data set containing 38 brain scans each with 20 total nodes.}

\usage{data(braingraphs)}

\format{
	The format is a data frame of 400 rows by 38 columns, with each column being a separate subject 
	and each row being a different edge between 2 nodes. Each column is a 20 by 20 matrix of brain 
	connections transformed into a vector. A value of 1 indicates that subject had a connection at that edge.
}

\keyword{datasets}