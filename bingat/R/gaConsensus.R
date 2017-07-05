gaConsensus <-
function(data, groups, iters=10, nRuns=1, popSize=200, method="manhattan", parallel=FALSE, cores=3){
	warning("This function has been deprecated.  Please use 'genAlg' or 'genAlgConsensus' instead.")
	ga <- genAlgConsensus(data, groups, .5, nRuns, parallel, cores, iters=iters, popSize=popSize, dataDist=method)
	return(ga)
}
