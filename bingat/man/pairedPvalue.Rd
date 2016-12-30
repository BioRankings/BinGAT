\name{pairedPvalue}
\alias{pairedPvalue}
\title{P-Value for Paired Data Results}

\description{This function returns the p-value of the significance of the difference in g-star values for paired data.}

\usage{pairedPvalue(data, type, groups, numPerms = 10, parallel = FALSE, cores = 3)}

\arguments{
  \item{data}{A data frame in which the columns (subjects) contain a 0/1 value for row (Node or Edge).}
  \item{type}{The type of graph being used (adjmatrix or adjmatrixlt).}
  \item{groups}{A vector of 0/1s used to denote the group membership of columns in the data.}
  \item{numPerms}{The number of permutations (ideally, 1000 or more).}
  \item{parallel}{TRUE or FALSE depending on whether the analysis will be parallelized for speed.}
  \item{cores}{The number of cores to use for parallelization. Ignored if parallel = FALSE.}
}

\value{A p-value.}

\author{Berkley Shands, Elena Deych, William D. Shannon}

\examples{
	data(braingraphs)
	
	grps <- c(rep(0, 19), rep(1, 19))
	numPermutations <- 1 ### This is set low for speed
	
	pval <- pairedPvalue(braingraphs, "adjMatrix", grps, numPermutations) 
	pval
}
