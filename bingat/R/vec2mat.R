vec2mat <-
function(x, type="adjMatrix"){
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
