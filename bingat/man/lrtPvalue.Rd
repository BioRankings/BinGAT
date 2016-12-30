\name{lrtPvalue}
\alias{lrtPvalue}
\title{Likelihood Ratio Test}

\description{This function returns the p-value of the significance between two groups.}

\usage{lrtPvalue(data, type, groups, numPerms = 10, parallel = FALSE, cores = 3)}

\arguments{
  \item{data}{A data frame in which the columns contain a 0/1 value for each node in the rows.}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{groups}{A vector of 0/1s used to denote the group membership of columns in the data.}
  \item{numPerms}{The number of permutations (ideally, 1000 or more).}
  \item{parallel}{TRUE or FALSE depending on whether the analysis will be parallelized for speed.}
  \item{cores}{The number of cores to use for parallelization. Ignored if parallel = FALSE.}
}

\value{The p-value for the difference between the two groups being tested.}

\author{Berkley Shands, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	grps <- sample(0:1, ncol(braingraphs), TRUE)
	numPermutations <- 1 ### This is set low for speed
	
	lrt <- lrtPvalue(braingraphs, "adjMatrix", grps, numPermutations) 
	lrt	
}
