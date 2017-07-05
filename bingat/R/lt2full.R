lt2full <-
function(x){	
	nodes <- getNumNodes(x, type="adjmatrixlt")
	
	y <- matrix(0, nodes, nodes)
	y[lower.tri(y)] <- x
	y <- y + t(y)
	
	return(as.vector(y))
}
