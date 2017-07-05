estTau <-
function(data, type, gstar){
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
