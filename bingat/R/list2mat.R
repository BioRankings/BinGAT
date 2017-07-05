list2mat <-
function(x, type="adjMatrix"){
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
