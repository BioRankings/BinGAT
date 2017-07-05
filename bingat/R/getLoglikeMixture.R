getLoglikeMixture <-
function(data, mixture, numConst){ 
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
