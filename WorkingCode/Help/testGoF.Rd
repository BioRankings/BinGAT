\name{testGoF}
\alias{testGoF}
\title{Test the Goodness of Fit}

\description{This function tests the goodness of fit for given a set of graphs.}

\usage{testGoF(data, type, numSims = 10, plot = TRUE, main)}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{numSims}{Number of simulations for Monte Carlo estimation of p-value(ideally, 1000 or more). Ignored if Chi-Square method is used.}
  \item{plot}{A boolean to create a plot of the results or not.}
  \item{main}{A title for the plot.}
}

\value{
	A list containing information about the goodness of fit and potentially a plot.  The list contains
	the Pearson statistics, degrees of freedom, and p-value, the G statistics and p-value, the Chi Squared
	statistics and p-value and finally the table with the observed and expected counts.
}

\examples{
	data(braingraphs)
	
	numSims <- 1 ### This is set low for speed
	braingof <- testGoF(braingraphs, "adjMatrix", numSims)
}
