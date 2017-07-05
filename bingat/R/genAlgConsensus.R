genAlgConsensus <-
function(data, covars, consensus=.5, numRuns=10, parallel=FALSE, cores=3, ...){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	
	if(consensus <= 0 || consensus > 1)
		stop("consensus must be greater than 0 and equal or less than 1")
	
	# Run the GA X times
	if(parallel){
		cl <- parallel::makeCluster(min(cores, numRuns)) 
		doParallel::registerDoParallel(cl)
		
		tryCatch({
					gaRes <- foreach::foreach(i=1:numRuns, .combine=list, .multicombine=TRUE, .inorder=FALSE, .packages=c("vegan", "HMP")) %dopar%{
						tempResults <- genAlg(data, covars, plot=FALSE, verbose=FALSE, ...)
						return(tempResults)
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		gaRes <- vector("list", numRuns)
		for(i in 1:numRuns)
			gaRes[[i]] <- genAlg(data, covars, plot=FALSE, verbose=FALSE, ...)
	}
	
	# Get all the best solutions
	bestSols <- sapply(gaRes, function(x){x$solutions[1,]})
	
	# Get the consensus solution vector
	consSol <- (rowSums(bestSols) >= (numRuns * consensus)) * 1
	
	# Get the selected Index's
	selInd <- which(consSol == 1)
	
	return(list(solutions=bestSols, consSol=consSol, selectedIndex=selInd))
}
