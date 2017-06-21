\name{plotMDS}
\alias{plotMDS}
\title{Plot MDS}

\description{This function plots all the data on an MDS plot.}

\usage{
	plotMDS(dataList, groups, estGstar = TRUE, paired = FALSE, returnCoords = FALSE, ..., data=NULL)
}

\arguments{
  \item{dataList}{A list where each element is a data frame in which the columns (subjects) 
  		contain a 0/1 value for row (Node or Edge).}
  \item{groups}{Deprecated. Each data set should be an element in dataList.}
  \item{estGstar}{When TRUE, the g star for every group is calculated and plotted.}
  \item{paired}{When TRUE, line segments between pairs will be drawn.}
  \item{returnCoords}{When TRUE, the MDS x-y coordinates will be returned.}
  \item{...}{Arguments to be passed to the plot method.}
  \item{data}{Deprecated. Replaced with dataList for clarity.}
}

\value{An MDS plot and if returnCoords is TRUE, a 2 column data frame containing the x-y coordinates of the data points is also returned.}

\author{Berkley Shands, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	### Break our data into two groups
	dataList <- list(braingraphs[,1:19], braingraphs[,20:38])
	
	### Basic plot
	plotMDS(dataList, main="MDS Plot")
	
	### Paired Plot
	plotMDS(dataList, paired=TRUE, main="Paired MDS Plot")
}
