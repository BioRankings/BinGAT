getGibbsMixture <-
function(data, type, desiredGroups, maxIter=50, digits=3){ 
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
