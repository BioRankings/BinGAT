gaCreation <-
function(data, popSize){
	ZERO_TO_ONE_RATIO <- 10 # Ratio of 0 to 1s for the random data
	SUGGESTION_COUNT <- 10 # Number starting points we should make from the data
	
	size <- ncol(data)
	population <- matrix(NA, popSize, size)
	
	# Make 10 starting points as long as our popSize is > 10
	if(popSize >= SUGGESTION_COUNT){
		# Get a rough starting point
		rstart <- apply(data, 2, mean)
		
		# Use the rough difference to make starting solutions
		breaks <- seq(.05, 1, 1/SUGGESTION_COUNT)
		suggestions <- matrix(0, length(breaks), length(rstart))
		for(i in 1:length(breaks))
			suggestions[i,] <- ifelse(rstart >= stats::quantile(rstart, breaks[i]), 1, 0)
		
		population[1:SUGGESTION_COUNT,] <- suggestions
		numCreated <- SUGGESTION_COUNT
	}else{
		numCreated <- 0
	}
	
	# Fill any remaining population spots with random solutions
	if(popSize != SUGGESTION_COUNT){
		for(child in (numCreated+1):popSize) 
			population[child,] <- sample(c(rep(0, ZERO_TO_ONE_RATIO), 1), size, replace=TRUE)
	}
	
	return(population)
}
