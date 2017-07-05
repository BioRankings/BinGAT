genAlg <-
function(data, covars, iters=50, popSize=200, earlyStop=0, dataDist="manhattan", covarDist="gower", 
		verbose=FALSE, plot=TRUE, minSolLen=NULL, maxSolLen=NULL){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	
	# Check for any bad numbers
	if(iters <= 0)
		stop("iters must be an integer greater than 0")
	if(popSize <= 0)
		stop("popSize must be an integer greater than 0")
	if(earlyStop < 0)
		stop("earlyStop must be an integer greater than or equal to 0")
	
	# Check distances
	if(dataDist != "manhattan")
		stop("data.dist must be manhattan.")
	if(covarDist != "euclidean" && covarDist != "gower")
		stop("covars.dist must be euclidean or gower.")
	
	# Define size
	size <- nrow(data)
	
	# Not ready for use yet
	penalty <- FALSE
	
	# Check stopping rules
	if(!is.null(minSolLen))
		if(minSolLen < 0 || minSolLen >= size)
			stop("minSolLen must be 0 or greater and less than the number of columns in data.")
	if(!is.null(maxSolLen))
		if(maxSolLen <= 0 || maxSolLen > size)
			stop("maxSolLen must be greater than 0 and less than or equal to the number columns in data.")
	if(!is.null(maxSolLen) && !is.null(minSolLen))
		if(maxSolLen < minSolLen)
			stop("maxSolLen must be bigger than minSolLen.")
	
	# Rotate the data
	data <- t(data)
	
	# Define some variables for use in the GA loop
	mutationChance <- 1/(size+1)
	elitism <- floor(popSize/5)
	evalSumm <- matrix(NA, iters, 6)
	newPopSize <- popSize - elitism
	newPopulation <- matrix(NA, newPopSize, size)
	parentProb <- stats::dnorm(1:popSize, mean=0, sd=(popSize/3))
	
	if(verbose){
		print("X. Current Step : Current Time Taken")
		runningTime <- proc.time()
		print(paste("1. Calculating Distances:", round((proc.time() - runningTime)[3], 3)))
	}
	
	# Set up our base distance matrix
	covarDists <- vegan::vegdist(covars, covarDist)
	
	# Get each columns distance contribution
	colDists <- vector("list", ncol(data))
	for(i in 1:ncol(data))
		colDists[[i]] <- vegan::vegdist(data[,i], dataDist)
	
	if(verbose)
		print(paste("2. Creating Starting Data:", round((proc.time() - runningTime)[3], 3)))
	
	# Create our starting data
	population <- gaCreation(data, popSize)
	
	if(verbose)
		print(paste("3. Scoring Starting Data:", round((proc.time() - runningTime)[3], 3)))
	
	# Score and sort
	evalVals <- rep(NA, popSize)
	for(e in 1:popSize)
		evalVals[e] <- gaScoring(population[e,], covarDists, colDists, dataDist, penalty, minSolLen, maxSolLen)
	population <- population[order(evalVals, decreasing=TRUE),]
	bestScoreValue <- max(evalVals)
	bestScoreCounter <- 0
	
	if(verbose)
		print(paste("4. Running Iterations:", round((proc.time() - runningTime)[3], 3)))
	
	# Run GA
	ptr <- proc.time()
	for(i in 1:iters){
		if(verbose){
			if(i %% round(iters/10) == 0)
				print(paste("Iteration - ", i, ": ", round((proc.time() - runningTime)[3], 3), sep=""))
		}
		# Cross over to fill rest of new population
		for(child in 1:newPopSize){
			parentIDs <- sample(1:popSize, 2, prob=parentProb)
			parents <- population[parentIDs,]
			crossOverPoint <- sample(0:size, 1)
			if(crossOverPoint == 0){
				newPopulation[child,] <- parents[2,]
			}else if(crossOverPoint == size){
				newPopulation[child,] <- parents[1,]
			}else{
				newPopulation[child,] <- c(parents[1,][1:crossOverPoint], parents[2,][(crossOverPoint+1):size])
			}
		}
		
		# Mutate all but elite
		if(mutationChance > 0){
			population[(elitism+1):popSize,] <- apply(newPopulation, 2, function(x){ifelse(stats::runif(newPopSize) < mutationChance, 1-x, x)})
		}else{
			population[(elitism+1):popSize,] <- newPopulation
		}
		
		# Score and sort our new solutions
		for(e in 1:popSize)
			evalVals[e] <- gaScoring(population[e,], covarDists, colDists, dataDist, penalty, minSolLen, maxSolLen)
		population <- population[order(evalVals, decreasing=TRUE),]
		evalSumm[i,] <- summary(evalVals)
		
		# Check if we want to stop early
		if(bestScoreValue == max(evalVals)){
			bestScoreCounter <- bestScoreCounter + 1
		}else{
			bestScoreCounter <- 0
			bestScoreValue <- max(evalVals)
		}
		
		if(bestScoreCounter == earlyStop && earlyStop != 0)
			break
	}	
	gaTime <- (proc.time() - ptr)[3]
	
	if(verbose)
		print(paste("5. Prettying Results", round((proc.time() - runningTime)[3], 3)))
	
	# Pretty up our data for returning
	rownames(population) <- paste("Solution", 1:nrow(population))
	colnames(population) <- colnames(data)
	rownames(evalSumm) <- paste("Iteration", 1:nrow(evalSumm))
	colnames(evalSumm) <- c("Best", "25%ile", "Median", "Mean", "75%ile", "Worst")
	
	evalVals <- matrix(evalVals[order(evalVals, decreasing=TRUE)], 1, length(evalVals))
	colnames(evalVals) <- paste("Solution ", 1:length(evalVals))
	rownames(evalVals) <- "Score"
	
	# Get selected columns using a consensus
	selIndex <- which(population[1,] == 1)
	sel <- colnames(data)[selIndex]
	
	# Get the nonselected columns
	nonSel <- colnames(data)[-selIndex]
	
	# Plot scoring summary
	if(plot)
		gaPlot(evalSumm)
	
	return(list(scoreSumm=evalSumm, solutions=population, scores=evalVals, time=gaTime, selected=sel, nonSelected=nonSel, selectedIndex=selIndex))
}
