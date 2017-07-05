full2lt <-
function(x){
	nodes <- getNumNodes(x, type="adjmatrix")
	
	y <- matrix(x, nodes, nodes, byrow=FALSE)
	y <- y[lower.tri(y)] 
	
	return(as.vector(y))
}
