gaScoring <-
function(indices, covarDists, colDists, distType, penalty, minSolLen, maxSolLen) {
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
