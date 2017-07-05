array2mat <-
function(x, type="adjMatrix"){	
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
