\name{plotMDS}
\alias{plotMDS}
\title{Plot MDS}

\description{This function plots all the data on an MDS plot.}

\usage{
	plotMDS(data, groups, estGstar = TRUE, paired = FALSE, returnCoords = FALSE, ...)
}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{groups}{A vector of integers used to denote the group membership of columns in the data.}
  \item{estGstar}{When TRUE, the g star for every group is calculated and plotted.}
  \item{paired}{When TRUE, line segments between pairs will be drawn.}
  \item{returnCoords}{When TRUE, the MDS x-y coordinates will be returned.}
  \item{...}{Arguments to be passed to the plot method.}
}

\value{An MDS plot and if returnCoords is TRUE, a 2 column data frame containing the x-y coordinates of the data points is also returned.}

\author{Berkley Shands, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	grps <- c(rep(0, 19), rep(1, 19))
	
	### Basic plot
	plotMDS(braingraphs, grps, main="My MDS Plot")
	
	### Paired Plot
	plotMDS(braingraphs, grps, TRUE, main="My Paired MDS Plot")
}
