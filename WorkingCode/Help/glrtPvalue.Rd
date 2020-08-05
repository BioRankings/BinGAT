\name{glrtPvalue}
\alias{glrtPvalue}
\title{GLRT Regression Results}

\description{This function returns the p-value of the significance of b1 in the regression model.}

\usage{
	glrtPvalue(dataList, type, groups, numPerms = 10, parallel = FALSE, 
		cores = 3, data = NULL)
}

\arguments{
  \item{dataList}{A list where each element is a data frame in which the columns (subjects) 
  		contain a 0/1 value for row (Node or Edge).}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{groups}{Deprecated. Each data set should be an element in dataList.}
  \item{numPerms}{Number of permutations.  In practice this should be at least 1,000.}
  \item{parallel}{TRUE or FALSE depending on whether the analysis will be parallelized for speed.}
  \item{cores}{The number of cores to use for parallelization. Ignored if parallel = FALSE.}
  \item{data}{Deprecated. Replaced with dataList for clarity.}
}

\value{A list containing the p-value, tau, logliklihood value, glrt value, bob1, b0, b1 and the hamming errors.}

\examples{
	data(braingraphs)
	
	### Break our data into two groups
	dataList <- list(braingraphs[,1:19], braingraphs[,20:38])
	
	### We use 1 for speed, should be at least 1,000
	numPerms <- 1
	
	res <- glrtPvalue(dataList, "adjMatrix", numPerms=numPerms) 
	res$pvalue
}
