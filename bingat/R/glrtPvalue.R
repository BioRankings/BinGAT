glrtPvalue <-
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
		
		# Turn the data sets into a list
		data1 <- dataList[,groups == unique(groups)[1], drop=FALSE]
		data2 <- dataList[,groups == unique(groups)[2], drop=FALSE]
		dataList <- list(data1, data2)
	}
	
	dataComb <- cbind(dataList[[1]], dataList[[2]])
	numSubTot <- ncol(dataComb)
	numSub1 <- ncol(dataList[[1]])
	
	# Run the glrt on the real data
	groups <- c(rep(0, numSub1), rep(1, ncol(dataList[[2]])))
	glrtObs <- glrtReg(dataComb, type, groups)
	
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		tryCatch({
					glrtPermute <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
						# Permute group membership
						samps <- sample(numSubTot, numSub1)
						permGrps <- rep(0, numSubTot)
						permGrps[samps] <- 1
						
						# Run glrt on permuted data
						glrt <- glrtReg(dataComb, type, permGrps)
						return(glrt)
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		glrtPermute <- rep(0, numPerms)
		for(i in 1:numPerms){
			# Permute group membership
			samps <- sample(numSubTot, numSub1)
			permGrps <- rep(0, numSubTot)
			permGrps[samps] <- 1
			
			# Run glrt on permuted data
			glrtPermute[i] <- glrtReg(dataComb, type, permGrps)
		}
	}
	
	# Calculate pvalue
	pvalue <- (sum(glrtPermute >= glrtObs) + 1)/(numPerms + 1)
	
	results <- c(glmReg(dataComb, type, groups), GLRT=glrtObs, pvalue=pvalue)
	return(results)
}
