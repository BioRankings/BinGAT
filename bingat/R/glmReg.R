glmReg <-
function(data, type, groups){	
	if(length(groups) == 1){
		# Get the gstar
		gstar <- estGStar(data)
		
		b0 <- gstar
		b0b1 <- NULL
		b1 <- NULL
		hammingError <- NULL
	}else{	
		# Estimate the gstars for each group
		b0 <- estGStar(data[,groups == 0, drop=FALSE])  
		b0b1 <- estGStar(data[,groups == 1, drop=FALSE])  
		b1 <- xor(b0, b0b1) * 1
		
		# estimate hamming error
		hammingError <- rep(0, ncol(data))
		for(i in 1:ncol(data))
			hammingError[i] <- calcDistance(data[,i], xor(b0, b1*groups[i]) * 1, type)
		
		# Get the gstars
		gstar <- matrix(0, nrow(data), ncol(data))
		for(i in 1:length(groups))
			gstar[,i] <- xor(b0, b1*groups[i])
	}
	
	# Get the tau and loglik values
	tau <- estTau(data, type, gstar)
	loglik <- estLogLik(data, type, gstar, tau)
	
	results <- list(b0, b1, b0b1, hammingError, loglik, tau)
	names(results) <- c("b0.covs0", "b1.Differences", "b0b1.covs1", "hammingError", "loglik", "tau")
	return(results)
}
