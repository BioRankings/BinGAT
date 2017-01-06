\name{glrtPvalue}
\alias{glrtPvalue}
\title{GLRT Regression Results}

\description{This function returns the p-value of the significance of b1 in the regression model.}

\usage{glrtPvalue(data, type, groups, numPerms = 10, parallel = FALSE, cores = 3)}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{groups}{A vector of 0/1s used to denote the group membership of columns in the data.}
  \item{numPerms}{The number of permutations (ideally, 1000 or more).}
  \item{parallel}{TRUE or FALSE depending on whether the analysis will be parallelized for speed.}
  \item{cores}{The number of cores to use for parallelization. Ignored if parallel = FALSE.}
}

\value{A list containing the results from \code{glm.reg} and \code{GLRT.reg} as well as the p-value.}

\author{Terrence Brooks, Berkley Shands, Skye Buckner-Petty, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	grps <- sample(0:1, ncol(braingraphs), TRUE)
	numPermutations <- 1 ### This is set low for speed
	
	glrt <- glrtPvalue(braingraphs, "adjMatrix", grps, numPermutations) 
	glrt
}
