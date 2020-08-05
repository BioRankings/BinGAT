

#The 'type' of data can refer to either:
#adjMatrix (m): 	An adjancency matrix made into a single vector
#adjMatrixLT (lt):	The upper or lower triangle of an adjancency matrix made into a single vector
#diag (d):  		The diagonal vector of an adj matrix 


library(network) 		# graphnetworkplot
library(matrixStats) 	# mixture models - logSumExp
library(gplots) 		# heatmap plotting
library(doParallel)		# paralleling
library(vegan)			# ga calcs


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### External
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ~~~~~~~~~~~~~~~~~~~~~
### pvalue functions
### ~~~~~~~~~~~~~~~~~~~~~
glrtPvalue <- function(dataList, type, groups, numPerms=10, parallel=FALSE, cores=3, data=NULL){
	# Check if data is still being used
	if(!is.null(data)){
		warning("'data' is deprecated. It has been replaced with dataList. View the help files for details.")
		dataList <- data
	}
	
	if(missing(dataList) || missing(type))
		stop("dataList and/or type is missing.")
	
	if(numPerms <= 0)
		stop("The number of permutations must be an integer greater than 0.")
	
	# Check if data/dataList isn't a list
	if(class(dataList) != "list"){
		warning("'dataList' should be a list of length 2. View the help files for details.")
		if(missing(groups))
			stop("groups is missing.")
		
		if(length(unique(groups)) != 2)
			stop("There must be exactly two groups.")
		
		# Turn the data sets into a list
		data1 <- dataList[,groups == unique(groups)[1], drop=FALSE]
		data2 <- dataList[,groups == unique(groups)[2], drop=FALSE]
		dataList <- list(data1, data2)
	}
	
	dataComb <- cbind(dataList[[1]], dataList[[2]])
	numSubTot <- ncol(dataComb)
	numSub1 <- ncol(dataList[[1]])
	
	# Run the glrt on the real data
	groups <- c(rep(0, numSub1), rep(1, ncol(dataList[[2]])))
	glrtObs <- glrtReg(dataComb, type, groups)
	
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		tryCatch({
					glrtPermute <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
						# Permute group membership
						samps <- sample(numSubTot, numSub1)
						permGrps <- rep(0, numSubTot)
						permGrps[samps] <- 1
						
						# Run glrt on permuted data
						glrt <- glrtReg(dataComb, type, permGrps)
						return(glrt)
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		glrtPermute <- rep(0, numPerms)
		for(i in 1:numPerms){
			# Permute group membership
			samps <- sample(numSubTot, numSub1)
			permGrps <- rep(0, numSubTot)
			permGrps[samps] <- 1
			
			# Run glrt on permuted data
			glrtPermute[i] <- glrtReg(dataComb, type, permGrps)
		}
	}
	
	# Calculate pvalue
	pvalue <- (sum(glrtPermute >= glrtObs) + 1)/(numPerms + 1)
	
	results <- c(glmReg(dataComb, type, groups), GLRT=glrtObs, pvalue=pvalue)
	return(results)
}

pairedPvalue  <- function(dataList, type, groups, numPerms=10, parallel=FALSE, cores=3, data=NULL){	
	# Check if data is still being used
	if(!is.null(data)){
		warning("'data' is deprecated. It has been replaced with dataList. View the help files for details.")
		dataList <- data
	}
	
	if(missing(dataList) || missing(type))
		stop("dataList and/or type is missing.")
	
	if(numPerms <= 0)
		stop("The number of permutations must be an integer greater than 0.")
	
	# Check if data/dataList isn't a list
	if(class(dataList) != "list"){
		warning("'dataList' should be a list of length 2. View the help files for details.")
		if(missing(groups))
			stop("groups is missing.")
		
		if(length(unique(groups)) != 2)
			stop("There must be exactly two groups.")
		
		if(max(groups) != 1 || min(groups) != 0)
			stop("'groups' must use 0 and 1 to denote groups.")
		
		data1 <- dataList[,groups == unique(groups)[1], drop=FALSE]
		data2 <- dataList[,groups == unique(groups)[2], drop=FALSE]
		
		dataList <- list(data1, data2)
	}
	
	numSub <- ncol(dataList[[1]])
	if(numSub != ncol(dataList[[2]]))
		stop("Both data sets must be the same size.")
	
	dataComb <- cbind(dataList[[1]], dataList[[2]])
	
	# Get our base distance
	gstarDistance <- calcDistance(estGStar(dataList[[1]]), estGStar(dataList[[2]]), type) 
	
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		tryCatch({
					permDistances <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
						# Permute the group membership
						samps <- sample(0:1, numSub, replace=TRUE)
						samps <- samps*numSub + 1:numSub
						
						gstar1 <- estGStar(dataComb[,samps, drop=FALSE])
						gstar2 <- estGStar(dataComb[,-samps, drop=FALSE])
						
						newDist <- calcDistance(gstar1, gstar2, type)
						return(newDist)
					}	
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		permDistances <- rep(0, numPerms)
		for(i in 1:numPerms){ 
			# Permute the group membership
			samps <- sample(0:1, numSub, replace=TRUE)
			samps <- samps*numSub + 1:numSub
			
			gstar1 <- estGStar(dataComb[,samps, drop=FALSE])
			gstar2 <- estGStar(dataComb[,-samps, drop=FALSE])
			
			permDistances[i] <- calcDistance(gstar1, gstar2, type)
		}
	}
	
	# Calculate pvalue
	pval <- (sum(permDistances >= gstarDistance) + 1)/(numPerms + 1) 
	return(pval)
}

lrtPvalue <- function(dataList, type, groups, numPerms=10, parallel=FALSE, cores=3, data=NULL){
	# Check if data is still being used
	if(!is.null(data)){
		warning("'data' is deprecated. It has been replaced with dataList. View the help files for details.")
		dataList <- data
	}
	
	if(missing(dataList) || missing(type))
		stop("dataList and/or type is missing.")
	
	if(numPerms <= 0)
		stop("The number of permutations must be an integer greater than 0.")
	
	# Check if data/dataList isn't a list
	if(class(dataList) != "list"){
		warning("'dataList' should be a list of length 2. View the help files for details.")
		if(missing(groups))
			stop("groups is missing.")
		
		if(length(unique(groups)) != 2)
			stop("There must be exactly two groups.")
		
		if(max(groups) != 1 || min(groups) != 0)
			stop("'groups' must use 0 and 1 to denote groups.")
		
		data1 <- dataList[,groups == unique(groups)[1], drop=FALSE]
		data2 <- dataList[,groups == unique(groups)[2], drop=FALSE]
		
		dataList <- list(data1, data2)
	}
	
	dataComb <- cbind(dataList[[1]], dataList[[2]])
	numSubTot <- ncol(dataComb)
	numSub1 <- ncol(dataList[[1]])
	
	# Get our logliks for each data set
	llC <- estLogLik(dataComb, type)
	ll1 <- estLogLik(dataList[[1]], type)
	ll2 <- estLogLik(dataList[[2]], type)
	
	# Get the lrt for the real data
	lrtObs <- -2 * (llC - ll1 - ll2)
	
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		tryCatch({
					lrtPermute  <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
						# Permute the samples
						samps <- sample(numSubTot, numSub1)
						
						# Get new groups log liks
						ll1 <- estLogLik(dataComb[,samps, drop=FALSE], type)
						ll2 <- estLogLik(dataComb[,-samps, drop=FALSE], type)
						
						lrt <- -2 * (llC - ll1 - ll2)
						return(lrt)
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		lrtPermute <- rep(0, numPerms)
		for(i in 1:numPerms){
			# Permute the samples
			samps <- sample(numSubTot, numSub1)
			
			# Get new groups log liks
			ll1 <- estLogLik(dataComb[,samps, drop=FALSE], type)
			ll2 <- estLogLik(dataComb[,-samps, drop=FALSE], type)
			
			lrtPermute[i] <- -2 * (llC - ll1 - ll2)
		}
	}
	
	# Calculate pvalue
	pValue <- (sum(lrtPermute > lrtObs) + 1)/(numPerms + 1)
	return(pValue)
}



### ~~~~~~~~~~~~~~~~~~~~~
### calc info functions
### ~~~~~~~~~~~~~~~~~~~~~
calcDistance <- function(x, y, type="", method="hamming"){
	if(missing(x) || missing(y))
		stop("x and/or y is missing.")
	
	if(tolower(method) == "hamming"){
		ret <- sum(as.logical(unlist(x)) != as.logical(unlist(y)))
		
		if(tolower(type) == "adjmatrix" || tolower(type) == "m") #Divide by 2 because of duplicates
			ret <- ret / 2
	}else{
		stop(sprintf("%s is unknown.", as.character(method)))
	}
	
	return(ret)
}

getNumNodes <- function(data, type){	
	if(missing(data) || missing(type))
		stop("data and/or type is missing.")
	
	data <- as.matrix(data)
	
	if(tolower(type) == "adjmatrix" || tolower(type) == "m"){
		nodes <- sqrt(nrow(data))
	}else if(tolower(type) == "adjmatrixlt" || tolower(type) == "lt"){
		c <- nrow(data)
		nodes <- (1 + sqrt(1 + 8*c))/2
	}else if(tolower(type) == "diag" || tolower(type) == "d"){
		nodes <- nrow(data)
	}else{
		stop(sprintf("%s is unknown.", as.character(type)))
	}
	
	return(nodes)	
}

getNumEdges <- function(nodes, type){
	if(missing(nodes) || missing(type))
		stop("nodes and/or type is missing.")
	
	if(tolower(type) == "diag" || tolower(type) == "d"){
		edges <- nodes
	}else if(tolower(type) == "adjmatrix" || tolower(type) == "adjmatrixlt"  || tolower(type) == "m" || tolower(type) == "lt"){
		edges <- choose(nodes, 2)
	}else{
		stop(sprintf("%s is unknown.", as.character(type)))
	}
	
	return(edges)
}

estGStar <- function(data){
	if(missing(data))
		stop("data is missing.")
	
	mean <- rowSums(data)/ncol(data)
	gstar <- 1*(mean > .5)
	
	return(gstar)
}

estTau <- function(data, type, gstar){
	if(missing(data) || missing(type) || missing(gstar))
		stop("data, type, and/or gstar is missing.")
	
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	
	if(!is.null(dim(gstar))){ # Check gstar is a single vector or matrix
		distToGStar <- NULL
		for(i in 1:ncol(gstar)){
			distToGStar[i] <- calcDistance(data[,i], gstar[,i], type)}
	}else{
		distToGStar <- apply(data, 2, function(x, g, type) {
					calcDistance(x, g, type)
				}, g=gstar, type=type)
	}
	
	sumDist <- sum(distToGStar)
	sumDist <- ifelse(sumDist==0, .Machine$double.xmin, sumDist) # Adjust sumDist if it is 0
	
	# Calculate tau
	num <- sumDist/(ncol(data) * edges)
	den <- 1-num
	tau <- -log(num/den)
	
	# Adjust for a negative tau in rare cases
	tau <- ifelse(tau < 0, 0, tau)
	
	return(tau)
}

estLogLik <- function(data, type, gstar, tau, g=NULL){
	if(missing(data) || missing(type))
		stop("data and/or type  is missing.")
	
	# Check if g is still used
	if(!is.null(g)){
		warning("'g' is deprecated. It has been replaced with gstar. View the help files for details.")
		gstar <- g
	}
	
	# Decide if we need to calculate gstar and tau too
	if(missing(gstar) && missing(tau)){ # Missing both tau and gstar
		mle <- estMLE(data, type)
		gstar <- mle$gstar
		tau <- mle$tau
	}else if(missing(gstar) || missing(tau)){ # Missing of of tau and gstar
		stop("tau or gstar is missing.")
	}
	
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	normConst <- 1 + exp(-tau)    
	normConst <- ifelse(normConst==0, .Machine$double.xmin, normConst) # Adjust normConst if it is 0
	
	# Calculate the distance from every point to the gstar
	if(!is.null(dim(gstar))){ # Check gstar is a single vector or matrix
		distToGStar <- rep(0, ncol(gstar))
		for(i in 1:ncol(gstar)){
			distToGStar[i] <- calcDistance(data[,i], gstar[,i], type)}
	}else{
		distToGStar <- apply(data, 2, function(x, g, type) {
					calcDistance(x, g, type)
				}, g=gstar, type=type)
	}
	
	# Calculate the log lik
	logLik <- -edges * ncol(data) * log(normConst) - tau * sum(distToGStar)
	return(logLik)
}

estMLE <- function(data, type){
	gstar <- estGStar(data)
	tau <- estTau(data, type, gstar)
	
	return(list(gstar=gstar, tau=tau))
}



### ~~~~~~~~~~~~~~~~~~~~~
### create data functions
### ~~~~~~~~~~~~~~~~~~~~~
rGibbs <- function(gstar, tau, type, numGraphs=1){
	if(missing(gstar) || missing(type) || missing(tau))
		stop("gstar, tau, prob, and/or type is missing.")
	
	if(numGraphs <= 0)
		stop("numGraphs must be an integer greater than 0.")
	
	nodes <- getNumNodes(gstar, type)
	edges <- getNumEdges(nodes, type)
	gstar <- as.vector(as.matrix(gstar))
	
	# Transform to a LT if we have a full matrix
	if(tolower(type) == "adjmatrix" || tolower(type) == "m")
		gstar <- full2lt(gstar)
	
	# Calculate the prob for a graph k-distance away and then sample distances from that
	possEdges <- 0:edges
	distProb <- exp(lgamma(edges+1) - lgamma(edges-possEdges+1) - lgamma(possEdges+1) - tau*possEdges-edges * log(1+exp(-tau)))
	gdists <- sample(possEdges, numGraphs, replace=TRUE, prob=distProb)
	
	# Take our starting graph and randomly flip-flop k nodes to get a graph k-distance away
	genData <- matrix(gstar, length(gstar), numGraphs)
	for(i in 1:numGraphs){ 
		nodesChange <- sample(1:length(gstar), gdists[i], replace=FALSE)
		genData[nodesChange, i] <- 1 * !genData[nodesChange, i]
	}
	
	# Transform back to a full matrix if we started with one
	if(tolower(type) == "adjmatrix" || tolower(type) == "m")
		genData <- apply(genData, 2, function(x){lt2full(x)})
	
	return(as.data.frame(genData)) 
}



### ~~~~~~~~~~~~~~~~~~~~~
### test data functions
### ~~~~~~~~~~~~~~~~~~~~~
testGoF <- function(data, type, numSims=10, plot=TRUE,  main){
	if(missing(data) || missing(type))
		stop("data and/or type is missing.")
	
	if(numSims <= 0)
		stop("The number of simulations must be an integer greater than 0.")
	
	numGraphs <- ncol(data)
	gstar <- estGStar(data)
	tau <- estTau(data, type, gstar)
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	
	# Calculate the prob for a graph k-distance away and calculate the expected counts of trees k-distance from gstar
	possEdges <- 0:edges
	distProb <- exp(lgamma(edges+1) - lgamma(edges-possEdges+1) - lgamma(possEdges+1) - tau*possEdges-edges * log(1+exp(-tau)))
	expCounts <- cbind(possEdges, numGraphs*distProb)
	
	# Calculate the observed counts of trees k-distance from gstar
	distToGStar <- apply(data, 2, function(x, g, type){calcDistance(x, g, type)}, g=gstar, type=type)
	distTable <- table(distToGStar)
	obsCounts <- cbind(as.integer(names(distTable)), as.integer(distTable))
	
	# Combine observed and expected counts into a single object
	obsExpCounts <- merge(expCounts, obsCounts, by=1, all=TRUE)
	colnames(obsExpCounts) <- c("dist", "expected", "observed")
	obsExpCounts[is.na(obsExpCounts)] <- 0
	methodA <- "Chisq"
	
	# Combine distances with a theoretical count < 5 (1 if there aren't enough at 5)
	signifExpect <- obsExpCounts$expected >= 5 
	if(sum(signifExpect) <= 1){
		signifExpect <- obsExpCounts$expected >= 1
		if(sum(signifExpect) < 1){
			stop("Expected counts below the threshold of 1. Goodness of fit cannot be calculated")
		}else{
			methodA <- "MC"
			warning("Expected counts below the threshold of 5; threshold lowered to 1 \n 
							P-value computed using Monte-Carlo simulation instead of asymptotic distribution.")
		}
	}
	
	# Determine main concentration of graphs
	sumSignifExp <- sum(signifExpect)
	lowerBound <- which(cumsum(signifExpect)==1)
	upperBound <- length(signifExpect)-which(cumsum(rev(signifExpect))==1)+1
	signifExpect[c(lowerBound, upperBound)] <- FALSE 
	
	# Collapse groups below/above the bounds
	lowerCount <- colSums(obsExpCounts[1:lowerBound, c("expected", "observed")])
	upperCount <- colSums(obsExpCounts[upperBound:length(signifExpect), c("expected", "observed")])
	mainCount <- obsExpCounts[signifExpect, c("expected", "observed")]
	
	# Bind all groups into 1 table
	oTable <- rbind(lowerCount, mainCount, upperCount)
	rownames(oTable) <- c(paste("<=", obsExpCounts[lowerBound, 1]), obsExpCounts[signifExpect, 1], paste(">=", obsExpCounts[upperBound, 1]))
	
	# Compute Pearson Chi-Squared Statistics
	pearsonStats <- sum(((oTable$observed - oTable$expected)^2)/oTable$expected)
	df <- nrow(oTable) - 1
	
	#Compute the Chi-Squared Statistics
	chisq <- chisq.test(oTable$observed, p=oTable$expected/numGraphs, simulate.p.value=ifelse(methodA == "MC", TRUE, FALSE), B=numSims)
	pvalueC <- chisq$p.value
	
	# Compute the G Statistics if needed
	if(sum(oTable$observed == 0) == 0){
		gStats <- 2*sum(oTable$observed * log(oTable$observed/oTable$expected))
		gdf <- nrow(oTable) - 1
		if(gdf <= 0){
			gdf <- nrow(oTable)
			warning("df is zero or negative; df replaced by the number of cells. (Conservative Test)")
		}
		pvalueG <- pchisq(gStats, df=gdf, ncp=0, lower.tail=FALSE, log.p=FALSE)
	}else{
		gStats <- NA
		pvalueG <- NA
	}
	
	pval <- ifelse(is.na(pvalueG), pvalueC, pvalueG)
	ptype <- ifelse(is.na(pvalueG), ifelse(methodA=="MC", "Monte-Carlo simulation", "Pearson Chi-square"), "G-test Statistics")
	
	# Plot the observed vs expected data
	if(plot){
		mycolor <- c("red", "blue")
		mylegend <- c("Expected", "Observed")
		signifExpect[c(lowerBound, upperBound)] <- TRUE 
		dist <- obsExpCounts[signifExpect, 1]
		matplot(dist, oTable, pch=19, type="p", col=mycolor)
		legend("topright", legend=mylegend, col=mycolor, pch=19)
		if(missing(main))
			title(c(paste(ptype), paste("P-value:", round(pval, 2))))
		else
			title(main)
	}
	
	results	<- list(ptype, df, pval, oTable)
	names(results) <- c("Method", "df", "pvalue",  "table")
	return(results)
}

getGibbsMixture <- function(data, type, desiredGroups, maxIter=50, digits=3){ 
	if(missing(data) || missing(desiredGroups) || missing(type))
		stop("data, type, and/or desiredGroups is missing.")
	
	if(maxIter <= 0)
		stop("maxIter must be an integer greater than 0.")
	
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	numSub <- ncol(data)
	
	# Set our starting points for gstars/weights/taus
	gstarsNew <- matrix(0, nrow(data), desiredGroups)
	weightsNew <- rep(0, desiredGroups)
	tausNew <- rep(0, desiredGroups)
	
	# Set our starting points for several checks
	gcheck <- FALSE
	tcheck <- FALSE
	wcheck <- FALSE
	converge <- FALSE
	
	# Randomly assign a starting group
	groups <- sample(desiredGroups, numSub, replace=TRUE)
	# Find the starting gstars/weights/taus for the given groups
	for(j in 1:desiredGroups){
		gstarsNew[,j] <- estGStar(data[, groups==j, drop=FALSE]) 
		weightsNew[j] <- sum(groups==j) / length(groups)
		tausNew[j] <- estTau(data[, groups==j, drop=FALSE], type, gstarsNew[,j])
	}
	
	for(iter in 1:maxIter){
		# Copy over our gstars/weights/taus
		gstarsOld <- gstarsNew
		weightsOld <- weightsNew
		tausOld <- tausNew
		
		# Calculate distances to the gstars
		distances <- matrix(NA, numSub, desiredGroups)
		for(j in 1:desiredGroups)
			distances[,j] <- apply(data, 2, function(x){calcDistance(x, gstarsNew[,j], type)})
		
		t <- matrix(NA, numSub, desiredGroups)
		pij <- matrix(NA, numSub, desiredGroups)
		for(i in 1:numSub){
			for(j in 1:desiredGroups)
				t[i, j] <- -edges*log(1 + exp(-tausNew[j])) - tausNew[j]*distances[i, j] + log(weightsNew[j])
			for(j in 1:desiredGroups)
				pij[i, j] <- 1/sum(exp(t[i,] - t[i, j]))
		}	
		
		# Check for covergence and quit if true
		if(gcheck && tcheck && wcheck){
			for(i in 1:numSub)
				groups[i] <- which(pij[i,] == max(pij[i,]))
			
			converge <- TRUE
			break
		}
		
		# Recompute gstars/weights/taus for new groups
		for(j in 1:desiredGroups){
			gnum <- apply(t(data) * pij[,j], 2, sum)
			gdem <- sum(pij[,j])
			gstarsNew[,j] <- ifelse(gnum/gdem > .5, 1, 0)  
			distj <- apply(data, 2, function(x){calcDistance(x, gstarsNew[,j], type)})
			
			weightsNew[j] <- sum(pij[,j])/numSub
			
			tnum <- sum(pij[,j] * distj)
			tdem <- sum(pij[,j] * edges) - tnum
			tfrac <- tnum/tdem
			tfrac <- ifelse(tfrac==0, .Machine$double.xmin, tfrac) # Adjust tfrac if it's 0
			tausNew[j] <- -log(tfrac) 
		}
		
		# Recalculate the checks for convergence
		wcheck <- all(round(weightsOld, digits) == round(weightsNew, digits))
		tcheck <- all(round(tausOld, digits) == round(tausNew, digits))
		if(wcheck && tcheck) # Only check gstars is tau and weight match
			gcheck <- all(gstarsOld == gstarsNew)
	}
	
	if(iter == maxIter && !converge)
		warning(sprintf("EM algorithm did not converge on %s groups with %s iterations", as.character(desiredGroups), as.character(maxIter)))
	
	return(list(weights=weightsNew, gstars=gstarsNew, taus=tausNew, converge=converge, 
					iterations=iter, numgroups=desiredGroups, type=type, pij=pij, group=groups))
}

getLoglikeMixture <- function(data, mixture, numConst){ 
	if(missing(data) || missing(mixture))
		stop("data and/or mixture is missing.")
	
	# Pull some variables from the mixture object
	type <- mixture$type	
	numGroups <- mixture$numgroups
	
	numSub <- ncol(data)
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	
	if(missing(numConst))
		numConst <- edges/nodes * 10 * log(numSub) # (2+edges) * log(numSub)
	
	# Calculate distances to the gstars
	distances <- matrix(NA, ncol(data), numGroups)
	for(j in 1:numGroups)
		distances[,j] <- apply(data, 2, function(x){calcDistance(x, mixture$gstars[[j]], type)})
	
	# Calculate inside of log sum exp equations
	temp <- matrix(NA, ncol(data), numGroups)
	for(i in 1:numSub)
		temp[i,] <- log(mixture$weights) - edges*log((1+exp(-mixture$taus))) - mixture$taus*distances[i,]
	
	# Calculate LL
	LL <- sum(apply(temp, 1, matrixStats::logSumExp))
	# Calculate BIC
	BIC <- -2*LL + numGroups * numConst
	
	return(list(ll=LL, bic=BIC))
}


### ~~~~~~~~~~~~~~~~~~~~~
### convert data functions
### ~~~~~~~~~~~~~~~~~~~~~
full2lt <- function(x){
	nodes <- getNumNodes(x, type="adjmatrix")
	
	y <- matrix(x, nodes, nodes, byrow=FALSE)
	y <- y[lower.tri(y)] 
	
	return(as.vector(y))
}

lt2full <- function(x){	
	nodes <- getNumNodes(x, type="adjmatrixlt")
	
	y <- matrix(0, nodes, nodes)
	y[lower.tri(y)] <- x
	y <- y + t(y)
	
	return(as.vector(y))
}

array2mat <- function(x, type="adjMatrix"){	
	nodes <- nrow(x)
	numSub <- dim(x)[3]
	
	if(tolower(type) == "adjmatrixlt" || tolower(type) == "lt"){
		y <- matrix(0, ((nodes^2 - nodes)/2), numSub)
		for(i in 1:numSub)
			y[,i] <- x[,,i][lower.tri(x[,,i])]
	}else{
		y <- matrix(0, nodes^2, numSub)
		for(i in 1:numSub)
			y[,i] <- x[,,i]
	}
	return(y)
}

list2mat <- function(x, type="adjMatrix"){
	nodes <- nrow(x[[1]])
	numSub <- length(x)
	
	if(tolower(type) == "adjmatrixlt" || tolower(type) == "lt"){
		y <- matrix(0, ((nodes^2 - nodes)/2), numSub)
		for(i in 1:numSub)
			y[,i] <- x[[i]][lower.tri(x[[i]])]
	}else{
		y <- matrix(0, nodes^2, numSub)
		for(i in 1:numSub)
			y[,i] <- as.vector(x[[i]])
	}
	return(y)
}

vec2mat <- function(x, type="adjMatrix"){
	nodes <- getNumNodes(x, type)
	
	if(tolower(type) == "adjmatrixlt" || tolower(type) == "lt"){
		y <- matrix(0, nodes, nodes)
		y[lower.tri(y)] <- x
		y <- y + t(y)
	}else if(tolower(type) == "diag" || tolower(type) == "d"){
		y <- matrix(0, nodes, nodes)
		diag(y) <- x
	}else{
		y <- matrix(x, nodes, nodes)
	}
	
	return(y)
}



### ~~~~~~~~~~~~~~~~~~~~~
### plotting functions
### ~~~~~~~~~~~~~~~~~~~~~
graphNetworkPlot <- function(data, type, main="Network Plot", labels, groupCounts, groupLabels){
	if(missing(data) || missing(type))
		stop("data and/or type is missing.")
	
	#Set up plot dot colors
	if(missing(groupCounts)){
		myColors <- "red"
	}else{
		allColors <- rainbow(length(groupCounts)+1)
		myColors <- NULL
		for(i in 1:length(groupCounts))
			myColors <- c(myColors, rep(allColors[i], groupCounts[i]))
	}
	
	# Take only the first column of data if it is multi columned
	if(class(data) == "data.frame" || class(data) == "matrix")
		data <- data[,1]
	
	dataVec <- vec2mat(data, type)
	graph <- network::network(as.matrix(dataVec), directed=FALSE)
	if(!missing(labels))
		network::network.vertex.names(graph) <- as.data.frame(labels)
	
	network::plot.network(graph, mode="circle", vertex.col=myColors, label=network::network.vertex.names(graph), main=main, edge.col="black")
	
	if(!missing(groupLabels))
		legend("topright", legend=groupLabels, fill=allColors, horiz=FALSE)
}

plotMDS <- function(dataList, groups, estGstar=TRUE, paired=FALSE, returnCoords=FALSE, ..., data=NULL){
	# Check if data is still being used
	if(!is.null(data)){
		warning("'data' is deprecated. It has been replaced with dataList. View the help files for details.")
		dataList <- data
	}
	
	if(missing(dataList))
		stop("dataList is missing.")
	
	# Check if data/dataList isn't a list
	if(class(dataList) != "list"){
		warning("'dataList' should be a list. View the help files for details.")
		if(missing(groups))
			stop("groups is missing.")
		
		# Turn the data sets into a list
		uniGroups <- unique(groups)
		dataTemp <- vector("list", length(uniGroups))
		for(i in 1:length(uniGroups))
			dataTemp[[i]] <- dataList[,groups == uniGroups[i], drop=FALSE]
		
		dataList <- dataTemp
	}
	
	# Get some variables
	numSubGrps <- sapply(dataList, ncol)
	numSub <- sum(numSubGrps)
	numGrps <- length(dataList)
	
	# Turn the data into a single data frame
	data <- do.call("cbind", dataList)
	
	# Make sure we have group names
	if(is.null(names(dataList))){
		grpNames <- paste("Data Set", 1:numGrps)
	}else{
		grpNames <- names(dataList)
	}
	
	# Check if we have paired data
	if(paired){
		if(numGrps != 2)
			stop("There must be exactly two data sets when paired = TRUE.")
		if(numSubGrps[1] != numSubGrps[2])
			stop("Both data sets must have the same number of subjects when paired = TRUE.")
	}
	
	# Set up the colors and sizes to use
	if(length(palette()) < numGrps)
		grDevices::palette(grDevices::rainbow(numGrps))
	myColors <- NULL
	for(i in 1:numGrps)
		myColors <- c(myColors, rep(i, numSubGrps[i]))
	myCexs <- rep(1, length(myColors))
	myPch <- rep(16, length(myColors))
	
	# Check if we need to calculate gstars
	if(estGstar){
		gstars <- matrix(0, nrow(dataList[[1]]), numGrps)
		for(i in 1:numGrps)
			gstars[,i] <- estGStar(dataList[[i]])
		colnames(gstars) <- grpNames
		
		data <- cbind(data, gstars)
		myColors <- c(myColors, grDevices::palette()[1:numGrps])
		myCexs <- c(myCexs, rep(1.5, numGrps))
		myPch <- c(myPch, rep(17, numGrps))
	}
	
	# Plot the mds
	loc <- cmdscale(dist(t(data), method="binary")) 
	plot(loc, xlab="MDS 1", ylab="MDS 2", pch=myPch, col=myColors, cex=myCexs, ...)
	
	# Add connecting lines if paired
	if(paired){
		groupSize <- numSubGrps[1]
		for(i in 1:groupSize)						
			segments(loc[i, 1], loc[i, 2], loc[i+groupSize, 1], loc[i+groupSize, 2], lty=3)
		
		if(estGstar)
			segments(loc[numSub+1, 1], loc[numSub+1, 2], loc[numSub+2, 1], loc[numSub+2, 2], lty=3) 	
	}
	
	# Label gstars
	if(estGstar)
		text(loc[(numSub+1):nrow(loc), 1], loc[(numSub+1):nrow(loc), 2], grpNames, pos=3)
	
	if(returnCoords)
		return(loc)
}

plotHeatmap <- function(data, type, names, ...){
	if(missing(data) || missing(type))
		stop("data and/or type is missing.")
	
	# Take only the first column of data if it is multi columned
	if(class(data) == "data.frame" || class(data) == "matrix")
		data <- data[,1]
	
	# Turn into a single matrix
	mat <- vec2mat(data, type)
	
	if(missing(names))
		names <- 1:ncol(mat)
	colnames(mat) <- names
	rownames(mat) <- names
	
	colfunc <- grDevices::colorRampPalette(c("blue", "grey"))
	
	gplots::heatmap.2(mat, symm=TRUE, Rowv=NA, dendrogram="none", trace="none", col=colfunc(10), ...)
}



### ~~~~~~~~~~~~~~~~~~~~~
### other functions
### ~~~~~~~~~~~~~~~~~~~~~
genAlg <- function(data, covars, iters=50, popSize=200, earlyStop=0, dataDist="manhattan", covarDist="gower", 
		verbose=FALSE, plot=TRUE, minSolLen=NULL, maxSolLen=NULL){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	
	# Check for any bad numbers
	if(iters <= 0)
		stop("iters must be an integer greater than 0")
	if(popSize <= 0)
		stop("popSize must be an integer greater than 0")
	if(earlyStop < 0)
		stop("earlyStop must be an integer greater than or equal to 0")
	
	# Check distances
	if(dataDist != "manhattan")
		stop("data.dist must be manhattan.")
	if(covarDist != "euclidean" && covarDist != "gower")
		stop("covars.dist must be euclidean or gower.")
	
	# Define size
	size <- nrow(data)
	
	# Not ready for use yet
	penalty <- FALSE
	
	# Check stopping rules
	if(!is.null(minSolLen))
		if(minSolLen < 0 || minSolLen >= size)
			stop("minSolLen must be 0 or greater and less than the number of columns in data.")
	if(!is.null(maxSolLen))
		if(maxSolLen <= 0 || maxSolLen > size)
			stop("maxSolLen must be greater than 0 and less than or equal to the number columns in data.")
	if(!is.null(maxSolLen) && !is.null(minSolLen))
		if(maxSolLen < minSolLen)
			stop("maxSolLen must be bigger than minSolLen.")
	
	# Rotate the data
	data <- t(data)
	
	# Define some variables for use in the GA loop
	mutationChance <- 1/(size+1)
	elitism <- floor(popSize/5)
	evalSumm <- matrix(NA, iters, 6)
	newPopSize <- popSize - elitism
	newPopulation <- matrix(NA, newPopSize, size)
	parentProb <- stats::dnorm(1:popSize, mean=0, sd=(popSize/3))
	
	if(verbose){
		print("X. Current Step : Current Time Taken")
		runningTime <- proc.time()
		print(paste("1. Calculating Distances:", round((proc.time() - runningTime)[3], 3)))
	}
	
	# Set up our base distance matrix
	covarDists <- vegan::vegdist(covars, covarDist)
	
	# Get each columns distance contribution
	colDists <- vector("list", ncol(data))
	for(i in 1:ncol(data))
		colDists[[i]] <- vegan::vegdist(data[,i], dataDist)
	
	if(verbose)
		print(paste("2. Creating Starting Data:", round((proc.time() - runningTime)[3], 3)))
	
	# Create our starting data
	population <- gaCreation(data, popSize)
	
	if(verbose)
		print(paste("3. Scoring Starting Data:", round((proc.time() - runningTime)[3], 3)))
	
	# Score and sort
	evalVals <- rep(NA, popSize)
	for(e in 1:popSize)
		evalVals[e] <- gaScoring(population[e,], covarDists, colDists, dataDist, penalty, minSolLen, maxSolLen)
	population <- population[order(evalVals, decreasing=TRUE),]
	bestScoreValue <- max(evalVals)
	bestScoreCounter <- 0
	
	if(verbose)
		print(paste("4. Running Iterations:", round((proc.time() - runningTime)[3], 3)))
	
	# Run GA
	ptr <- proc.time()
	for(i in 1:iters){
		if(verbose){
			if(i %% round(iters/10) == 0)
				print(paste("Iteration - ", i, ": ", round((proc.time() - runningTime)[3], 3), sep=""))
		}
		# Cross over to fill rest of new population
		for(child in 1:newPopSize){
			parentIDs <- sample(1:popSize, 2, prob=parentProb)
			parents <- population[parentIDs,]
			crossOverPoint <- sample(0:size, 1)
			if(crossOverPoint == 0){
				newPopulation[child,] <- parents[2,]
			}else if(crossOverPoint == size){
				newPopulation[child,] <- parents[1,]
			}else{
				newPopulation[child,] <- c(parents[1,][1:crossOverPoint], parents[2,][(crossOverPoint+1):size])
			}
		}
		
		# Mutate all but elite
		if(mutationChance > 0){
			population[(elitism+1):popSize,] <- apply(newPopulation, 2, function(x){ifelse(stats::runif(newPopSize) < mutationChance, 1-x, x)})
		}else{
			population[(elitism+1):popSize,] <- newPopulation
		}
		
		# Score and sort our new solutions
		for(e in 1:popSize)
			evalVals[e] <- gaScoring(population[e,], covarDists, colDists, dataDist, penalty, minSolLen, maxSolLen)
		population <- population[order(evalVals, decreasing=TRUE),]
		evalSumm[i,] <- summary(evalVals)
		
		# Check if we want to stop early
		if(bestScoreValue == max(evalVals)){
			bestScoreCounter <- bestScoreCounter + 1
		}else{
			bestScoreCounter <- 0
			bestScoreValue <- max(evalVals)
		}
		
		if(bestScoreCounter == earlyStop && earlyStop != 0)
			break
	}	
	gaTime <- (proc.time() - ptr)[3]
	
	if(verbose)
		print(paste("5. Prettying Results", round((proc.time() - runningTime)[3], 3)))
	
	# Pretty up our data for returning
	rownames(population) <- paste("Solution", 1:nrow(population))
	colnames(population) <- colnames(data)
	rownames(evalSumm) <- paste("Iteration", 1:nrow(evalSumm))
	colnames(evalSumm) <- c("Best", "25%ile", "Median", "Mean", "75%ile", "Worst")
	
	evalVals <- matrix(evalVals[order(evalVals, decreasing=TRUE)], 1, length(evalVals))
	colnames(evalVals) <- paste("Solution ", 1:length(evalVals))
	rownames(evalVals) <- "Score"
	
	# Get selected columns using a consensus
	selIndex <- which(population[1,] == 1)
	sel <- colnames(data)[selIndex]
	
	# Get the nonselected columns
	nonSel <- colnames(data)[-selIndex]
	
	# Plot scoring summary
	if(plot)
		gaPlot(evalSumm)
	
	return(list(scoreSumm=evalSumm, solutions=population, scores=evalVals, time=gaTime, selected=sel, nonSelected=nonSel, selectedIndex=selIndex))
}

genAlgConsensus <- function(data, covars, consensus=.5, numRuns=10, parallel=FALSE, cores=3, ...){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	
	if(consensus <= 0 || consensus > 1)
		stop("consensus must be greater than 0 and equal or less than 1")
	
	# Run the GA X times
	if(parallel){
		cl <- parallel::makeCluster(min(cores, numRuns)) 
		doParallel::registerDoParallel(cl)
		
		tryCatch({
					gaRes <- foreach::foreach(i=1:numRuns, .combine=list, .multicombine=TRUE, .inorder=FALSE, .packages=c("vegan", "HMP")) %dopar%{
						tempResults <- genAlg(data, covars, plot=FALSE, verbose=FALSE, ...)
						return(tempResults)
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		gaRes <- vector("list", numRuns)
		for(i in 1:numRuns)
			gaRes[[i]] <- genAlg(data, covars, plot=FALSE, verbose=FALSE, ...)
	}
	
	# Get all the best solutions
	bestSols <- sapply(gaRes, function(x){x$solutions[1,]})
	
	# Get the consensus solution vector
	consSol <- (rowSums(bestSols) >= (numRuns * consensus)) * 1
	
	# Get the selected Index's
	selInd <- which(consSol == 1)
	
	return(list(solutions=bestSols, consSol=consSol, selectedIndex=selInd))
}



### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Internal
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ~~~~~~~~~~~~~~~~~~~~~
### pval help functions
### ~~~~~~~~~~~~~~~~~~~~~
glmReg <- function(data, type, groups){	
	if(length(groups) == 1){
		# Get the gstar
		gstar <- estGStar(data)
		
		b0 <- gstar
		b0b1 <- NULL
		b1 <- NULL
		hammingError <- NULL
	}else{	
		# Estimate the gstars for each group
		b0 <- estGStar(data[,groups == 0, drop=FALSE])  
		b0b1 <- estGStar(data[,groups == 1, drop=FALSE])  
		b1 <- xor(b0, b0b1) * 1
		
		# estimate hamming error
		hammingError <- rep(0, ncol(data))
		for(i in 1:ncol(data))
			hammingError[i] <- calcDistance(data[,i], xor(b0, b1*groups[i]) * 1, type)
		
		# Get the gstars
		gstar <- matrix(0, nrow(data), ncol(data))
		for(i in 1:length(groups))
			gstar[,i] <- xor(b0, b1*groups[i])
	}
	
	# Get the tau and loglik values
	tau <- estTau(data, type, gstar)
	loglik <- estLogLik(data, type, gstar, tau)
	
	results <- list(b0, b1, b0b1, hammingError, loglik, tau)
	names(results) <- c("b0.covs0", "b1.Differences", "b0b1.covs1", "hammingError", "loglik", "tau")
	return(results)
}

glrtReg <- function(data, type, groups){
	regCoeffHA <- glmReg(data, type, groups)
	regCoeffH0 <- glmReg(data, type, 0)
	
	glrt <- 2 * (regCoeffHA$loglik - regCoeffH0$loglik)
	return(glrt)
}



### ~~~~~~~~~~~~~~~~~~~~~
### ga functions
### ~~~~~~~~~~~~~~~~~~~~~
gaConsensus <- function(data, groups, iters=10, nRuns=1, popSize=200, method="manhattan", parallel=FALSE, cores=3){
	warning("This function has been deprecated.  Please use 'genAlg' or 'genAlgConsensus' instead.")
	ga <- genAlgConsensus(data, groups, .5, nRuns, parallel, cores, iters=iters, popSize=popSize, dataDist=method)
	return(ga)
}

gaScoring <- function(indices, covarDists, colDists, distType, penalty, minSolLen, maxSolLen) {
	BAD_RETURN <- -2 # Return worse than cor could do
	
	numSel <- sum(indices)
	
	# Check if nothing is selected
	if(numSel == 0) 
		return(BAD_RETURN) 
	# Check if we dont have enough selected
	if(!is.null(minSolLen)) 
		if(numSel < minSolLen)
			return(BAD_RETURN)
	# Check if we dont have too many selected
	if(!is.null(maxSolLen)) 
		if(numSel > maxSolLen)
			return(BAD_RETURN) 
	
	edges <- which(indices==1)
	combinedSolDists <- Reduce("+", colDists[edges])
	
	# Get the correlation and penalize it based on the number of columns selected
	mycor <- stats::cor(combinedSolDists, covarDists)
	if(penalty)
		mycor <- mycor * (length(indices)-sum(indices))/(length(indices)-1)
	
	return(mycor)
}

gaCreation <- function(data, popSize){
	ZERO_TO_ONE_RATIO <- 10 # Ratio of 0 to 1s for the random data
	SUGGESTION_COUNT <- 10 # Number starting points we should make from the data
	
	size <- ncol(data)
	population <- matrix(NA, popSize, size)
	
	# Make 10 starting points as long as our popSize is > 10
	if(popSize >= SUGGESTION_COUNT){
		# Get a rough starting point
		rstart <- apply(data, 2, mean)
		
		# Use the rough difference to make starting solutions
		breaks <- seq(.05, 1, 1/SUGGESTION_COUNT)
		suggestions <- matrix(0, length(breaks), length(rstart))
		for(i in 1:length(breaks))
			suggestions[i,] <- ifelse(rstart >= stats::quantile(rstart, breaks[i]), 1, 0)
		
		population[1:SUGGESTION_COUNT,] <- suggestions
		numCreated <- SUGGESTION_COUNT
	}else{
		numCreated <- 0
	}
	
	# Fill any remaining population spots with random solutions
	if(popSize != SUGGESTION_COUNT){
		for(child in (numCreated+1):popSize) 
			population[child,] <- sample(c(rep(0, ZERO_TO_ONE_RATIO), 1), size, replace=TRUE)
	}
	
	return(population)
}

gaPlot <- function(evalSumm){
	graphics::plot(evalSumm[,4], type="l", ylab="Score", ylim=c(0, 1), lwd=2, main="Eval Scores by Iteration", xlab="Iteration")
	graphics::lines(evalSumm[,6], col="red", lwd=2)
	graphics::lines(evalSumm[,1], col="blue", lwd=2)
	graphics::legend("topleft", colnames(evalSumm)[c(4, 6, 1)], pch=16, col=c("black", "red", "blue"))
}



### ~~~~~~~~~~~~~~~~~~~~~
### other/new
### ~~~~~~~~~~~~~~~~~~~~~
getGibbsMixtureNew <- function(data, type, desiredGroups, maxIter=50, digits=3){ 
	if(missing(data) || missing(desiredGroups) || missing(type))
		stop("data, type, and/or desiredGroups is missing.")
	
	if(maxIter <= 0)
		stop("maxIter must be an integer greater than 0.")
	
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	numSub <- ncol(data)
	
	# Set our starting points for gstars/weights/taus
	gstarsNew <- matrix(0, nrow(data), desiredGroups)
	weightsNew <- rep(0, desiredGroups)
	tausNew <- rep(0, desiredGroups)
	
	# Set our starting points for several checks
	gcheck <- FALSE
	tcheck <- FALSE
	wcheck <- FALSE
	converge <- FALSE
	
	# Randomly assign a starting group
	groups <- matrix(NA, numSub, maxIter+1)
	newGroups <- sample(desiredGroups, numSub, replace=TRUE)
	groups[,1] <- newGroups
	# Find the starting gstars/weights/taus for the given groups
	for(j in 1:desiredGroups){
		gstarsNew[,j] <- estGStar(data[, newGroups == j, drop=FALSE]) 
		weightsNew[j] <- sum(newGroups == j) / length(groups)
		tausNew[j] <- estTau(data[, newGroups == j, drop=FALSE], type, gstarsNew[,j])
	}
	
	for(iter in 1:maxIter){
		# Copy over our gstars/weights/taus
		gstarsOld <- gstarsNew
		weightsOld <- weightsNew
		tausOld <- tausNew
		
		# Calculate distances to the gstars
		distances <- matrix(NA, numSub, desiredGroups)
		for(j in 1:desiredGroups)
			distances[,j] <- apply(data, 2, function(x){calcDistance(x, gstarsNew[,j], type)})
		
		t <- matrix(NA, numSub, desiredGroups)
		pij <- matrix(NA, numSub, desiredGroups)
		for(i in 1:numSub){
			for(j in 1:desiredGroups)
				t[i, j] <- -edges*log(1 + exp(-tausNew[j])) - tausNew[j]*distances[i, j] + log(weightsNew[j])
			for(j in 1:desiredGroups)
				pij[i, j] <- 1/sum(exp(t[i,] - t[i, j]))
		}	
		
		newGroups <- apply(pij, 1, function(x) which(x == max(x)))
		groups[,iter+1] <- newGroups
		
		# Recompute gstars/weights/taus for new groups
		for(j in 1:desiredGroups){
			gstarsNew[,j] <- estGStar(data[, newGroups == j, drop=FALSE]) 
			weightsNew[j] <- sum(newGroups == j) / length(groups)
			tausNew[j] <- estTau(data[, newGroups == j, drop=FALSE], type, gstarsNew[,j])
		}
		
		# Check for covergence and quit if true
		if(all(newGroups == groups[,iter])){
			converge <- TRUE
			break
		}
	}
	
	if(iter == maxIter && !converge)
		warning(sprintf("EM algorithm did not converge on %s groups with %s iterations", as.character(desiredGroups), as.character(maxIter)))
	
	return(list(weights=weightsNew, gstars=gstarsNew, taus=tausNew, converge=converge, 
					iterations=iter, numgroups=desiredGroups, type=type, pij=pij, group=groups))
}




