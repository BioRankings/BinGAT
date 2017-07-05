lrtPvalue <-
function(dataList, type, groups, numPerms=10, parallel=FALSE, cores=3, data=NULL){
	# Check if data is still being used
	if(!is.null(data)){
		warning("'data' is deprecated. It has been replaced with dataList. View the help files for details.")
		dataList <- data
	}
	
	if(missing(dataList) || missing(type))
		stop("dataList and/or type is missing.")
	
	if(numPerms <= 0)
		stop("The number of permutations must be an integer greater than 0.")
	
	# Check if data/dataList isn't a list
	if(class(dataList) != "list"){
		warning("'dataList' should be a list of length 2. View the help files for details.")
		if(missing(groups))
			stop("groups is missing.")
		
		if(length(unique(groups)) != 2)
			stop("There must be exactly two groups.")
		
		if(max(groups) != 1 || min(groups) != 0)
			stop("'groups' must use 0 and 1 to denote groups.")
		
		data1 <- dataList[,groups == unique(groups)[1], drop=FALSE]
		data2 <- dataList[,groups == unique(groups)[2], drop=FALSE]
		
		dataList <- list(data1, data2)
	}
	
	dataComb <- cbind(dataList[[1]], dataList[[2]])
	numSubTot <- ncol(dataComb)
	numSub1 <- ncol(dataList[[1]])
	
	# Get our logliks for each data set
	llC <- estLogLik(dataComb, type)
	ll1 <- estLogLik(dataList[[1]], type)
	ll2 <- estLogLik(dataList[[2]], type)
	
	# Get the lrt for the real data
	lrtObs <- -2 * (llC - ll1 - ll2)
	
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		tryCatch({
					lrtPermute  <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
						# Permute the samples
						samps <- sample(numSubTot, numSub1)
						
						# Get new groups log liks
						ll1 <- estLogLik(dataComb[,samps, drop=FALSE], type)
						ll2 <- estLogLik(dataComb[,-samps, drop=FALSE], type)
						
						lrt <- -2 * (llC - ll1 - ll2)
						return(lrt)
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		lrtPermute <- rep(0, numPerms)
		for(i in 1:numPerms){
			# Permute the samples
			samps <- sample(numSubTot, numSub1)
			
			# Get new groups log liks
			ll1 <- estLogLik(dataComb[,samps, drop=FALSE], type)
			ll2 <- estLogLik(dataComb[,-samps, drop=FALSE], type)
			
			lrtPermute[i] <- -2 * (llC - ll1 - ll2)
		}
	}
	
	# Calculate pvalue
	pValue <- (sum(lrtPermute > lrtObs) + 1)/(numPerms + 1)
	return(pValue)
}
