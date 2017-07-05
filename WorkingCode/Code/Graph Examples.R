
source("Graph Functions.R")

braingraphs <- read.csv("braingraphs.csv", row.names=1)


### ~~~~~~~~~~~~~~~~~~~~~
### pvalue functions
### ~~~~~~~~~~~~~~~~~~~~~
glrtPvalue.Test <- function(){
	data(braingraphs)
	
	### Break our data into two groups
	dataList <- list(braingraphs[,1:19], braingraphs[,20:38])
	
	### We use 1 for speed, should be at least 1,000
	numPerms <- 1
	
	res <- glrtPvalue(dataList, "adjMatrix", numPerms=numPerms) 
	res$pvalue
}

pairedPvalue.Test <- function(){
	data(braingraphs)
	
	### Break our data into two groups
	dataList <- list(braingraphs[,1:19], braingraphs[,20:38])
	
	### We use 1 for speed, should be at least 1,000
	numPerms <- 1
	
	pval <- pairedPvalue(dataList, "adjMatrix", numPerms=numPerms) 
	pval
}

lrtPvalue.Test <- function(){
	data(braingraphs)
	
	### Break our data into two groups
	dataList <- list(braingraphs[,1:19], braingraphs[,20:38])
	
	### We use 1 for speed, should be at least 1,000
	numPerms <- 1
	
	lrt <- lrtPvalue(dataList, "adjMatrix", numPerms=numPerms) 
	lrt	
}



### ~~~~~~~~~~~~~~~~~~~~~
### calc info functions
### ~~~~~~~~~~~~~~~~~~~~~
calcDistance.Test <- function(){
	data(braingraphs)
	
	dist <- calcDistance(braingraphs[,1], braingraphs[,2], "adjMatrix")
	dist
}

getNumNodes.Test <- function(){
	data(braingraphs)
	
	brainnodes <- getNumNodes(braingraphs, "adjMatrix")
	brainnodes
}

getNumEdges.Test <- function(){
	data(braingraphs)
	
	brainnodes <- getNumNodes(braingraphs, "adjMatrix")
	brainedges <- getNumEdges(brainnodes, "adjMatrix")
	brainedges
}

estGStar.Test <- function(){
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	braingstar[1:25]
}

estTau.Test <- function(){
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	braintau
}

estLogLik.Test <- function(){
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	brainll <- estLogLik(braingraphs, "adjMatrix", braingstar, braintau)
	brainll
}

estMLE.Test <- function(){
	data(braingraphs)
	
	brainmle <- estMLE(braingraphs, "adjMatrix") 
	brainmle
}



### ~~~~~~~~~~~~~~~~~~~~~
### create data functions
### ~~~~~~~~~~~~~~~~~~~~~
rGibbs.Test <- function(){
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs)
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	randombraingraphs <- rGibbs(braingstar, braintau, "adjMatrix", 3) 
	randombraingraphs[1:5,]
}



### ~~~~~~~~~~~~~~~~~~~~~
### test data functions
### ~~~~~~~~~~~~~~~~~~~~~
testGoF.Test <- function(){
	data(braingraphs)
	
	numSims <- 1 ### This is set low for speed
	braingof <- testGoF(braingraphs, "adjMatrix", numSims)
}

getGibbsMixture.Test <- function(){
	data(braingraphs)
	
	braingm <- getGibbsMixture(braingraphs, "adjMatrix", 5) 
}

getLoglikeMixture.Test <- function(){
	data(braingraphs)
	
	braingm <- getGibbsMixture(braingraphs, "adjMatrix", 5)
	brainlm <- getLoglikeMixture(braingraphs, braingm)
	brainlm
	
	### By running the loglik mixture over several groups you can find which is the optimal
	\dontrun{
		mixtures <- NULL
		for(i in 1:5){
			tempgm <- getGibbsMixture(braingraphs, "adjMatrix", i)
			mixtures[i] <- getLoglikeMixture(braingraphs, tempgm)$bic
		}
		
		bestgroupnum <- which(min(mixtures) == mixtures)
		bestgroupnum
	}
}


### ~~~~~~~~~~~~~~~~~~~~~
### plotting functions
### ~~~~~~~~~~~~~~~~~~~~~
graphNetworkPlot.Test <- function(){
	data(braingraphs)
	
	main <- "Brain Connections"
	gc <- c(5, 5, 4, 6)
	gl <- c("Grp1", "Grp2", "Grp3", "Grp4")
	
	graphNetworkPlot(braingraphs[,1], "adjMatrix", main, groupCounts=gc, groupLabels=gl)
}

plotMDS.Test <- function(){
	data(braingraphs)
	
	### Break our data into two groups
	dataList <- list(braingraphs[,1:19], braingraphs[,20:38])
	
	### Basic plot
	plotMDS(dataList, main="MDS Plot")
	
	### Paired Plot
	plotMDS(dataList, paired=TRUE, main="Paired MDS Plot")
}

plotHeatmap.Test <- function(){
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	plotHeatmap(braingstar, "adjMatrix")	
}



### ~~~~~~~~~~~~~~~~~~~~~
### other functions
### ~~~~~~~~~~~~~~~~~~~~~
genAlgConsensus.Test <- function(){
	\dontrun{
		data(braingraphs)
		
		### Set covars to just be group membership
		covars <- matrix(c(rep(0, 19), rep(1, 19)))
		
		### We use low numbers for speed. The exact numbers to use depend
		### on the data being used, but generally the higher iters and popSize 
		### the longer it will take to run.  earlyStop is then used to stop the
		### run early if the results aren't improving.
		iters <- 500
		popSize <- 200
		earlyStop <- 250
		numRuns <- 3
		
		gaRes <- genAlgConsensus(braingraphs, covars, .5, numRuns, FALSE, 3, 
				iters, popSize, earlyStop)
	}
}

genAlg.Test <- function(){
	\dontrun{
		data(braingraphs)
		
		### Set covars to just be group membership
		covars <- matrix(c(rep(0, 19), rep(1, 19)))
		
		### We use low numbers for speed. The exact numbers to use depend
		### on the data being used, but generally the higher iters and popSize 
		### the longer it will take to run.  earlyStop is then used to stop the
		### run early if the results aren't improving.
		iters <- 500
		popSize <- 200
		earlyStop <- 250
		
		gaRes <- genAlg(braingraphs, covars, iters, popSize, earlyStop)
	}
}




