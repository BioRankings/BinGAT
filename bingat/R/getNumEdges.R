getNumEdges <-
function(nodes, type){
	if(missing(nodes) || missing(type))
		stop("nodes and/or type is missing.")
	
	if(tolower(type) == "diag" || tolower(type) == "d"){
		edges <- nodes
	}else if(tolower(type) == "adjmatrix" || tolower(type) == "adjmatrixlt"  || tolower(type) == "m" || tolower(type) == "lt"){
		edges <- choose(nodes, 2)
	}else{
		stop(sprintf("%s is unknown.", as.character(type)))
	}
	
	return(edges)
}
