plotMDS <-
function(dataList, groups, estGstar=TRUE, paired=FALSE, returnCoords=FALSE, ..., data=NULL){
	# Check if data is still being used
	if(!is.null(data)){
		warning("'data' is deprecated. It has been replaced with dataList. View the help files for details.")
		dataList <- data
	}
	
	if(missing(dataList))
		stop("dataList is missing.")
	
	# Check if data/dataList isn't a list
	if(class(dataList) != "list"){
		warning("'dataList' should be a list. View the help files for details.")
		if(missing(groups))
			stop("groups is missing.")
		
		# Turn the data sets into a list
		uniGroups <- unique(groups)
		dataTemp <- vector("list", length(uniGroups))
		for(i in 1:length(uniGroups))
			dataTemp[[i]] <- dataList[,groups == uniGroups[i], drop=FALSE]
		
		dataList <- dataTemp
	}
	
	# Get some variables
	numSubGrps <- sapply(dataList, ncol)
	numSub <- sum(numSubGrps)
	numGrps <- length(dataList)
	
	# Turn the data into a single data frame
	data <- do.call("cbind", dataList)
	
	# Make sure we have group names
	if(is.null(names(dataList))){
		grpNames <- paste("Data Set", 1:numGrps)
	}else{
		grpNames <- names(dataList)
	}
	
	# Check if we have paired data
	if(paired){
		if(numGrps != 2)
			stop("There must be exactly two data sets when paired = TRUE.")
		if(numSubGrps[1] != numSubGrps[2])
			stop("Both data sets must have the same number of subjects when paired = TRUE.")
	}
	
	# Set up the colors and sizes to use
	if(length(palette()) < numGrps)
		grDevices::palette(grDevices::rainbow(numGrps))
	myColors <- NULL
	for(i in 1:numGrps)
		myColors <- c(myColors, rep(i, numSubGrps[i]))
	myCexs <- rep(1, length(myColors))
	myPch <- rep(16, length(myColors))
	
	# Check if we need to calculate gstars
	if(estGstar){
		gstars <- matrix(0, nrow(dataList[[1]]), numGrps)
		for(i in 1:numGrps)
			gstars[,i] <- estGStar(dataList[[i]])
		colnames(gstars) <- grpNames
		
		data <- cbind(data, gstars)
		myColors <- c(myColors, grDevices::palette()[1:numGrps])
		myCexs <- c(myCexs, rep(1.5, numGrps))
		myPch <- c(myPch, rep(17, numGrps))
	}
	
	# Plot the mds
	loc <- cmdscale(dist(t(data), method="binary")) 
	plot(loc, xlab="MDS 1", ylab="MDS 2", pch=myPch, col=myColors, cex=myCexs, ...)
	
	# Add connecting lines if paired
	if(paired){
		groupSize <- numSubGrps[1]
		for(i in 1:groupSize)						
			segments(loc[i, 1], loc[i, 2], loc[i+groupSize, 1], loc[i+groupSize, 2], lty=3)
		
		if(estGstar)
			segments(loc[numSub+1, 1], loc[numSub+1, 2], loc[numSub+2, 1], loc[numSub+2, 2], lty=3) 	
	}
	
	# Label gstars
	if(estGstar)
		text(loc[(numSub+1):nrow(loc), 1], loc[(numSub+1):nrow(loc), 2], grpNames, pos=3)
	
	if(returnCoords)
		return(loc)
}
