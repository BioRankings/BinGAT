estLogLik <-
function(data, type, gstar, tau, g=NULL){
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
