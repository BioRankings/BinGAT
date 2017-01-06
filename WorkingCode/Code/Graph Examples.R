	
	
	
	source("Graph Functions.R")
	
	
	braingraphs <- read.csv("braingraphs.csv", row.names=1)

##########################################################
### getNumNodes
##########################################################
	data(braingraphs)
	
	brainnodes <- getNumNodes(braingraphs, "adjMatrix")
	brainnodes
	
##########################################################
### getNumEdges
##########################################################
	data(braingraphs)

	brainnodes <- getNumNodes(braingraphs, "adjMatrix")
	brainedges <- getNumEdges(brainnodes, "adjMatrix")
	brainedges
	
##########################################################
### estGStar
##########################################################
	data(braingraphs)

	braingstar <- estGStar(braingraphs) 
	braingstar[1:25]
	
##########################################################
### estTau
##########################################################
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	braintau
	
##########################################################
### estLogLik
##########################################################
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	brainll <- estLogLik(braingraphs, "adjMatrix", braingstar, braintau)
	brainll
	
##########################################################
### estMLE
##########################################################
	data(braingraphs)

	brainmle <- estMLE(braingraphs, "adjMatrix") 
	brainmle
	
##########################################################
### rGibbs
##########################################################
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs)
	braintau <- estTau(braingraphs, "adjMatrix", braingstar)
	randombraingraphs <- rGibbs(braingstar, braintau, "adjMatrix", 3) 
	randombraingraphs[1:5,]
	
##########################################################
### testGoF
##########################################################
	data(braingraphs)
	
	numSims <- 1 ### This is set low for speed
	braingof <- testGoF(braingraphs, "adjMatrix", numSims)
	
##########################################################
### getGibbsMixture
##########################################################
	data(braingraphs)

	braingm <- getGibbsMixture(braingraphs, "adjMatrix", 5) 
	
##########################################################
### getLoglikeMixture
##########################################################
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
	
##########################################################
### glrtPvalue
##########################################################
	data(braingraphs)
	
	grps <- sample(0:1, ncol(braingraphs), TRUE)
	numPermutations <- 1 ### This is set low for speed
	
	glrt <- glrtPvalue(braingraphs, "adjMatrix", grps, numPermutations) 
	glrt
	
##########################################################
### graphNetworkPlot
##########################################################
	data(braingraphs)
	
	main <- "Brain Connections"
	gc <- c(5, 5, 4, 6)
	gl <- c("Grp1", "Grp2", "Grp3", "Grp4")
	
	graphNetworkPlot(braingraphs[,1], "adjMatrix", main, groupCounts=gc, groupLabels=gl)
	
##########################################################
### calcDistance
##########################################################
	data(braingraphs)
	
	dist <- calcDistance(braingraphs[,1], braingraphs[,2], "adjMatrix")
	dist
	
##########################################################
### pairedPvalue
##########################################################
	data(braingraphs)
	
	grps <- c(rep(0, 19), rep(1, 19))
	numPermutations <- 1 ### This is set low for speed
	
	pval <- pairedPvalue(braingraphs, "adjMatrix", grps, numPermutations) 
	pval

##########################################################
### plotMDS
##########################################################
	data(braingraphs)
	
	grps <- c(rep(0, 19), rep(1, 19))
	
	### Basic plot
	plotMDS(braingraphs, grps, main="My MDS Plot")
	
	### Paired Plot
	plotMDS(braingraphs, grps, paired=TRUE, main="My Paired MDS Plot")
	
	### Contour Plot
	sicknessVal <- c(sample(5:10, 19, TRUE), sample(1:5, 19, TRUE))
	plotMDS(braingraphs, grps, FALSE, TRUE, sicknessVal, main="My Plot")
	
##########################################################
### plotHeatmap
##########################################################
	data(braingraphs)
	
	braingstar <- estGStar(braingraphs) 
	plotHeatmap(braingstar, "adjMatrix")	

##########################################################
### lrtPvalue
##########################################################
	data(braingraphs)
	
	grps <- sample(0:1, ncol(braingraphs), TRUE)
	numPermutations <- 1 ### This is set low for speed
	
	lrt <- lrtPvalue(braingraphs, "adjMatrix", grps, numPermutations) 
	lrt	
	
##########################################################
### ga
##########################################################
	data(braingraphs)
	
	grps <- c(rep(0, 19), rep(1, 19))
	iters <- 10 ### This is set low for speed
	nRuns <- 5 ### This is set low for speed
	
	consensus <- gaConsensus(braingraphs, grps, iters, nRuns) 
	consensus$corr[1:5]
	consensus$solutions[1:2,]

	
	
	
