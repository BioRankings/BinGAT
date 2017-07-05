pairedPvalue <-
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
	
	numSub <- ncol(dataList[[1]])
	if(numSub != ncol(dataList[[2]]))
		stop("Both data sets must be the same size.")
	
	dataComb <- cbind(dataList[[1]], dataList[[2]])
	
	# Get our base distance
	gstarDistance <- calcDistance(estGStar(dataList[[1]]), estGStar(dataList[[2]]), type) 
	
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		tryCatch({
					permDistances <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
						# Permute the group membership
						samps <- sample(0:1, numSub, replace=TRUE)
						samps <- samps*numSub + 1:numSub
						
						gstar1 <- estGStar(dataComb[,samps, drop=FALSE])
						gstar2 <- estGStar(dataComb[,-samps, drop=FALSE])
						
						newDist <- calcDistance(gstar1, gstar2, type)
						return(newDist)
					}	
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		permDistances <- rep(0, numPerms)
		for(i in 1:numPerms){ 
			# Permute the group membership
			samps <- sample(0:1, numSub, replace=TRUE)
			samps <- samps*numSub + 1:numSub
			
			gstar1 <- estGStar(dataComb[,samps, drop=FALSE])
			gstar2 <- estGStar(dataComb[,-samps, drop=FALSE])
			
			permDistances[i] <- calcDistance(gstar1, gstar2, type)
		}
	}
	
	# Calculate pvalue
	pval <- (sum(permDistances >= gstarDistance) + 1)/(numPerms + 1) 
	return(pval)
}
