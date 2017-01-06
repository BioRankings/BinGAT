

#The 'type' of data can refer to either:
#adjMatrix: 	An adjancency matrix made into a single vector
#adjMatrixLT:	The upper or lower triangle of an adjancency matrix made into a single vector
#diag:  		The diagonal vector of an adj matrix 


library(network) 		#for graphnetworkplot
library(matrixStats) 	#for mixture models - logSumExp
#library(vegan) 		#for contour plotting
#library(akima) 			#for contour plotting
library(gplots) 		#for heatmap plotting
library(genalg)			#for ga algorithim
library(doParallel)		#for paralleling

calcDistance <- function(x, y, type="", method="hamming"){
	if(missing(x) || missing(y))
		stop("x and/or y is missing.")
	
	if(tolower(method) == "hamming"){
		ret <- sum(as.logical(unlist(x)) != as.logical(unlist(y)))
		
		if(tolower(type) == "adjmatrix") #Divide by 2 because of duplicates
			ret <- ret / 2
	}else{
		stop(sprintf("%s is unknown.", as.character(method)))
	}
	
	return(ret)
}

getNumNodes <- function(data, type){	
	if(missing(data) || missing(type))
		stop("data and/or type is missing.")
	
	data <- as.matrix(data)
	
	if(tolower(type) == "adjmatrix"){
		nodes <- sqrt(nrow(data))
	}else if(tolower(type) == "adjmatrixlt"){
		c <- nrow(data)
		nodes <- (1+sqrt(1+8*c))/2
	}else if(tolower(type) == "diag"){
		nodes <- nrow(data)
	}else{
		stop(sprintf("%s is unknown.", as.character(type)))
	}
	
	return(nodes)	
}

getNumEdges <- function(nodes, type){
	if(missing(nodes) || missing(type))
		stop("nodes and/or type is missing.")
	
	if(tolower(type) == "diag"){
		edges <- nodes
	}else if(tolower(type) == "adjmatrix" || tolower(type) == "adjmatrixlt"){
		edges <- choose(nodes, 2)
	}else{
		stop(sprintf("%s is unknown.", as.character(type)))
	}
	
	return(edges)
}

estGStar <- function(data){
	if(missing(data))
		stop("data is missing.")
	
	mean <- rowSums(data)/ncol(data)
	gstar <- 1*(mean > .5)
	
	return(gstar)
}

estTau <- function(data, type, gstar){
	if(missing(data) || missing(type) || missing(gstar))
		stop("data, type, and/or gstar is missing.")
	
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	
	if(!is.null(dim(gstar))){ #Check gstar is a single vector or matrix
		distToGStar <- NULL
		for(i in 1:ncol(gstar)){
			distToGStar[i] <- calcDistance(data[,i], gstar[,i], type)}
	}else{
		distToGStar <- apply(data, 2, function(x, g, type) {calcDistance(x, g, type)}, g=gstar, type=type)
	}
	
	sumDist <- sum(distToGStar)
	sumDist <- ifelse(sumDist==0, .Machine$double.xmin, sumDist) #Adjust sumDist if it is 0
	
	num <- sumDist/(ncol(data) * edges)
	den <- 1-num
	tau <- -log(num/den)
	
	tau <- ifelse(tau < 0, 0, tau) #Adjust for a negative tau in rare cases
	
	return(tau)
}

estLogLik <- function(data, type, g, tau){
	if(missing(data) || missing(type) || missing(g) || missing(tau))
		stop("data, type, g, and/or tau is missing.")
	
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	normConst <- 1+exp(-tau)    
	###normConst <- ifelse(normConst==0, .Machine$double.xmin, normConst) #Adjust normConst if it is 0
	
	if(!is.null(dim(g))){ #Check gstar is a single vector or matrix
		distToGStar <- NULL
		for(i in 1:ncol(g)){
			distToGStar[i] <- calcDistance(data[,i], g[,i], type)}
	}else{
		distToGStar <- apply(data, 2, function(x, g, type) {calcDistance(x, g, type)}, g=g, type=type)
	}
	
	LogLik <- -edges * ncol(data) * log(normConst) - tau * sum(distToGStar)
	
	return(LogLik)
}

estMLE <- function(data, type){
	gstar <- estGStar(data)
	tau <- estTau(data, type, gstar)
	
	return(list(gstar=gstar, tau=tau))
}

rGibbs <- function(gstar, tau, type, numGraphs=1){
	if(missing(gstar) || missing(type) || missing(tau))
		stop("gstar, tau, prob, and/or type is missing.")
	
	if(numGraphs <= 0)
		stop("numGraphs must be an integer greater than 0.")
	
	nodes <- getNumNodes(gstar, type)
	edges <- getNumEdges(nodes, type)
	gstar <- as.vector(as.matrix(gstar))
	
	if(tolower(type) == "adjmatrix")
		gstar <- full2lt(gstar)
	
	#Calculate the prob for a graph k-distance away and then sample distances from that
	possEdges <- 0:edges
	distProb <- exp(lgamma(edges+1) - lgamma(edges-possEdges+1) - lgamma(possEdges+1) - tau*possEdges-edges * log(1+exp(-tau)))
	gdists <- sample(possEdges, numGraphs, replace=TRUE, prob=distProb)
	
	#Take our starting graph and randomly flip-flop k nodes to get a graph k-distance away
	genData <- matrix(gstar, length(gstar), numGraphs)
	for(i in 1:numGraphs){ 
		nodesChange <- sample(1:length(gstar), gdists[i], replace=FALSE)
		genData[nodesChange, i] <- 1*!genData[nodesChange, i]
	}
	
	if(tolower(type) == "adjmatrix")
		genData <- apply(genData, 2, function(x){lt2full(x)})
	
	return(as.data.frame(genData)) 
}
### updated to include just one p-value (asymptotic or permutations + title

testGoF <- function(data, type, numSims=10, plot=TRUE,  main){
	if(missing(data) || missing(type))
		stop("data and/or type is missing.")
	
	if(numSims <= 0)
		stop("The number of simulations must be an integer greater than 0.")
	
	numGraphs <- ncol(data)
	gstar <- estGStar(data)
	tau <- estTau(data, type, gstar)
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	
	#Calculate the prob for a graph k-distance away and calculate the expected counts of trees k-distance from gstar
	possEdges <- 0:edges
	distProb <- exp(lgamma(edges+1) - lgamma(edges-possEdges+1) - lgamma(possEdges+1) - tau*possEdges-edges * log(1+exp(-tau)))
	expCounts <- cbind(possEdges, numGraphs*distProb)
	
	#Calculate the observed counts of trees k-distance from gstar
	distToGStar <- apply(data, 2, function(x, g, type){calcDistance(x, g, type)}, g=as.vector(gstar), type=type)
	distTable <- table(distToGStar)
	obsCounts <- cbind(as.integer(names(distTable)), as.integer(distTable))
	
	#Combine observed and expected counts into a single object
	obsExpCounts <- merge(expCounts, obsCounts, by=1, all=TRUE)
	colnames(obsExpCounts) <- c("dist", "expected", "observed")
	obsExpCounts[is.na(obsExpCounts)] <- 0
	methodA <- "Chisq"
	
	#Combine distances with a theoretical count < 5 (1 if there aren't enough at 5)
	signifExpect <- obsExpCounts$expected >= 5 
	if(sum(signifExpect) <= 1){
		signifExpect <- obsExpCounts$expected >= 1
		if(sum(signifExpect) < 1){
			stop("Expected counts below the threshold of 1. Goodness of fit cannot be calculated")
		}else{
			methodA <- "MC"
			warning("Expected counts below the threshold of 5; threshold lowered to 1 \n 
							P-value computed using Monte-Carlo simulation instead of asymptotic distribution.")
		}
	}
	
	#Determine main concentration of graphs
	sumSignifExp <- sum(signifExpect)
	lowerBound <- which(cumsum(signifExpect)==1)
	upperBound <- length(signifExpect)-which(cumsum(rev(signifExpect))==1)+1
	signifExpect[c(lowerBound, upperBound)] <- FALSE 
	
	#Collapse groups below/above the bounds
	lowerCount <- colSums(obsExpCounts[1:lowerBound, c("expected", "observed")])
	upperCount <- colSums(obsExpCounts[upperBound:length(signifExpect), c("expected", "observed")])
	mainCount <- obsExpCounts[signifExpect, c("expected", "observed")]
	
	#Bind all groups into 1 table
	oTable <- rbind(lowerCount, mainCount, upperCount)
	rownames(oTable) <- c(paste("<=", obsExpCounts[lowerBound, 1]), obsExpCounts[signifExpect, 1], paste(">=", obsExpCounts[upperBound, 1]))
	
	#Compute Pearson Chi-Squared Statistics
	pearsonStats <- sum(((oTable$observed-oTable$expected)^2)/oTable$expected)
	df <- nrow(oTable)-1
	
	#Compute the Chi-Squared Statistics
	chisq <- chisq.test(oTable$observed, p=oTable$expected/numGraphs, simulate.p.value=ifelse(methodA=="MC",TRUE,FALSE), B=numSims)
	pvalueC <- chisq$p.value
	
	#Compute the G Statistics if needed
	if(sum(oTable$observed==0)==0){
		gStats <- 2*sum(oTable$observed * log(oTable$observed/oTable$expected))
		gdf <- nrow(oTable)-1
		if(gdf <= 0){
			gdf <- nrow(oTable)
			warning("df is zero or negative; df replaced by the number of cells. (Conservative Test)")
		}
		pvalueG <- pchisq(gStats, df=gdf, ncp=0, lower.tail=FALSE, log.p=FALSE)
	}else{
		gStats <- NA
		pvalueG <- NA
	}
	
	pval <- ifelse(is.na(pvalueG), pvalueC, pvalueG)
	ptype <- ifelse(is.na(pvalueG), ifelse(methodA=="MC", "Monte-Carlo simulation", "Pearson Chi-square"), "G-test Statistics")
	
	#Plot the observed vs expected data
	if(plot){
		mycolor <- c("red", "blue")
		mylegend <- c("Expected", "Observed")
		signifExpect[c(lowerBound, upperBound)] <- TRUE 
		dist <- obsExpCounts[signifExpect, 1]
		matplot(dist, oTable, pch=19, type="p", col=mycolor)
		legend("topright", legend=mylegend, col=mycolor, pch=19)
		if(missing(main))
			title(c(paste(ptype), paste("P-value:", round(pval, 2))))
		else
			title(main)
	}
	
	results	<- list(ptype, df, pval, oTable)
	names(results) <- c("Method", "df", "pvalue",  "table")
	return(results)
}

getGibbsMixture <- function(data, type, desiredGroups, maxIter=50, digits=3){ 
	if(missing(data) || missing(desiredGroups) || missing(type))
		stop("data, type, and/or desiredGroups is missing.")
	
	if(maxIter <= 0)
		stop("maxIter must be an integer greater than 0.")
	
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	groups <- sample(1:desiredGroups, ncol(data), replace=TRUE)
	
	gstars <- list()
	weights <- NULL
	taus <- NULL
	
	gcheck <- 0
	tcheck <- 0
	wcheck <- 1
	group <- NULL
	converge <- NA
	
	#Find the starting gstars/tau/weights for the given groups
	for(j in 1:desiredGroups){
		weights[j] <- sum(groups==j) / length(groups)
		gstars[[j]] <- estGStar(data[, groups==j, drop=FALSE]) 
		taus[j] <- estTau(data[, groups==j, drop=FALSE], type, gstars[[j]])
	}
	
	for(iter in 1:maxIter){
		#Calculate distances to the gstars
		distances <- matrix(NA, ncol(data), desiredGroups)
		for(j in 1:desiredGroups)
			distances[,j] <- apply(data, 2, function(x){calcDistance(x, gstars[[j]], type)})
		
		t <- matrix(NA, ncol(data), desiredGroups)
		pij <- matrix(NA, ncol(data), desiredGroups)
		for(i in 1:ncol(data)){
			for(j in 1:desiredGroups)
				t[i, j] <- -edges*log(1 + exp(-taus[j])) - taus[j]*distances[i, j] + log(weights[j])
			for(j in 1:desiredGroups)
				pij[i, j] <- 1/sum(exp(t[i,] - t[i, j]))
		}	
		
		#Check for covergence and quit if true
		if(gcheck + tcheck + wcheck == 0){
			for(i in 1:ncol(data))
				group <- c(group, which(pij[i,] == max(pij[i,])))
			
			converge <- TRUE
			break
		}
		
		#Recompute gstars/taus/weights for new groups
		gstars2 <- list()
		weights2 <- NULL
		taus2 <- NULL
		for(j in 1:desiredGroups){
			weights2[j] <- sum(pij[,j]) / ncol(data)
			gnum <- apply(t(data)*pij[,j], 2, sum)
			gdem <- sum(pij[,j])
			gstars2[[j]] <- ifelse(gnum/gdem > .5, 1, 0)  
			distj <- apply(data, 2, function(x){calcDistance(x, gstars2[[j]], type)})
			
			tnum <- sum(pij[,j] * distj)
			tdem <- sum(pij[,j] * edges) - tnum
			tfrac <- tnum/tdem
			tfrac <- ifelse(tfrac==0, .Machine$double.xmin, tfrac)
			
			taus2[j] <- -log(tfrac) 
		}
		
		#Recalculate the checks for convergence
		gcheck <- 0
		tcheck <- 0
		wcheck <- 0
		for(j in 1:desiredGroups){
			gcheck <- gcheck + sum(gstars[[j]] != gstars2[[j]])
			tcheck <- tcheck + as.numeric(round(taus[j], digits) != round(taus2[j], digits))
			wcheck <- wcheck + as.numeric(round(weights[j], digits) != round(weights2[j], digits))
		}
		
		if(iter == maxIter){
			warning(sprintf("EM algorithm did not converge on %s groups with %s iterations", as.character(desiredGroups), as.character(maxIter)))
			converge <- FALSE
		}
		
		gstars <- gstars2
		taus <- taus2
		weights <- weights2
	}
	
	return(list(weights=weights2, gstars=gstars2, taus=taus2, converge=converge, 
					iterations=iter, numgroups=desiredGroups, type=type, pij=pij, group=group))
}

getLoglikeMixture <- function(data, mixture, numConst){ 
	if(missing(data) || missing(mixture))
		stop("data and/or mixture is missing.")
	
	type <- mixture$type	
	numGroups <- mixture$numgroups
	numGraphs <- ncol(data)
	nodes <- getNumNodes(data, type)
	edges <- getNumEdges(nodes, type)
	
	#Calculate distances to the gstars
	distances <- matrix(NA, ncol(data), numGroups)
	for(j in 1:numGroups)
		distances[,j] <- apply(data, 2, function(x){calcDistance(x, mixture$gstars[[j]], type)})
	
	t <- matrix(NA, ncol(data), numGroups)
	for(i in 1:numGraphs)
		t[i,] <- log(mixture$weights) - edges*log((1+exp(-mixture$taus))) - mixture$taus*distances[i,]
	
	LL <- sum(apply(t, 1, matrixStats::logSumExp))
	
	if(missing(numConst))
		numConst <- edges/nodes * 10 * log(numGraphs) #(2+edges) * log(numGraphs)
	BIC <- -2*LL + numGroups * numConst
	
	return(list(ll=LL, bic=BIC))
}

# Converts full vectorized matrix into its vectorized lower triangle
full2lt <- function(x){
	nn <- getNumNodes(x, type="adjmatrix")
	y <- matrix(x, nn, nn, byrow=FALSE)
	y <- y[lower.tri(y)] 
	
	return(as.vector(y))
}

# Converts vectorized lower triangle of a matrix into its vectorized full matrix
lt2full <- function(x){	
	nn <- getNumNodes(x, type="adjmatrixlt")
	y <- matrix(0, nn, nn)
	y[lower.tri(y)] <- x
	y <- y + t(y)
	
	return(as.vector(y))
}

# Converts array of matrices to single matrices of adjacency matrix vectors
array2mat <- function(x, type="adjMatrix"){	
	if(tolower(type) == "adjmatrixlt"){
		y <- matrix(0, (((dim(x)[1]*dim(x)[1])-dim(x)[1])/2), dim(x)[3])
		for(i in 1:dim(x)[3])
			y[,i] <- x[,,i][lower.tri(x[,,i])]
	}else{
		y <- matrix(0, dim(x)[1]*dim(x)[1], dim(x)[3])
		for(i in 1:dim(x)[3])
			y[,i] <- x[,,i]
	}
	return(y)
}

# Converts lists of matrices to adjacency matrix
list2mat <- function(x, type="adjMatrix"){	
	if(tolower(type) == "adjmatrixlt"){
		y <- matrix(0, (((dim(x[[1]])[1]*dim(x[[1]])[1]) - dim(x[[1]])[1])/2), length(x))
		for(i in 1:length(x))
			y[,i] <- x[[i]][lower.tri(x[[i]])]
	}else{
		y <- matrix(0, dim(x[[1]])[1]*dim(x[[1]])[1], length(x))
		for(i in 1:length(x))
			y[,i] <- as.vector(x[[i]])
	}
	return(y)
}    

# Converts vectors of matrices (or lower/upper triangles of matrices to adjacency matrix
vec2mat <- function(x, type="adjMatrix"){
	y <- NULL
	nn <- getNumNodes(x, type)
	
	if(tolower(type) == "adjmatrixlt"){
		y <- matrix(0, nn, nn)
		y[lower.tri(y)] <- x
		y <- y + t(y)
	}else if(tolower(type) == "diag"){
		y <- matrix(0, nn, nn)
		diag(y) <- x
	}else{
		y <- matrix(x, nn, nn)
	}
	
	return(y)
}

glmReg <- function(data, type, groups){	
	if(dim(table(groups)) != 2){
		gstar <- estGStar(data)
		b0 <- gstar
		b0b1 <- NULL
		b1 <- NULL
		hammingError <- NULL
		tau <- estTau(data, type, gstar)	
		loglik <- estLogLik(data, type, gstar, tau)	
	}else{	
		#Estimate the gstars for each group
		b0 <- estGStar(data[,groups==unique(groups)[1]])  
		b0b1 <- estGStar(data[,groups==unique(groups)[2]])  
		b1 <- xor(b0, b0b1)*1
		
		index <- as.matrix(1:length(groups))
		hammingError <- apply(index, 1, function(x, b0, b1, grps, datap, typ){
					calcDistance(datap[,x], xor(b0, b1*grps[x]), typ)
				}, b0=b0, b1=b1, grps=groups, datap=data, typ=type)
		
		#Get the tau and loglik values
		gstar <- matrix(0, nrow(data), ncol(data))
		for(i in 1:length(groups))
			gstar[,i] <- as.matrix(xor(b0, b1*groups[i])) 
		
		tau <- estTau(data, type, gstar)
		loglik <- estLogLik(data, type, gstar, tau)
	}
	
	results <- list(b0, b1, b0b1, hammingError, loglik, tau)
	names(results) <- c("b0.covs0", "b1.Differences", "b0b1.covs1", "hammingError", "loglik", "tau")
	
	return(results)
}

glrtReg <- function(data, type, groups){	
	regCoeffHA <- glmReg(data, type, groups)
	loglikHA <- regCoeffHA$loglik
	
	regCoeffH0 <- glmReg(data, type, rep(0, ncol(data)))
	loglikH0 <- regCoeffH0$loglik
	
	glrt <- 2 * (loglikHA-loglikH0)
	return(glrt)
}

glrtPvalue <- function(data, type, groups, numPerms=10, parallel=FALSE, cores=3){
	if(missing(data) || missing(type) || missing(groups))
		stop("data, type, and/or groups is missing.")
	
	if(numPerms <= 0)
		stop("The number of permutations must be an integer greater than 0.")
	
	if(length(unique(groups)) != 2)
		stop("There must be exactly two groups.")
	
	if(max(groups) != 1 || min(groups) != 0)
		stop("'groups' must use 0 and 1 to denote groups.")
	
	x1 <- sum(groups==unique(groups)[1])
	x2 <- sum(groups==unique(groups)[2])
	
	if(x1 <= 0 || x2 <= 0)
		stop("Each group must have at least 1 subject.")
	
	reg <- glmReg(data, type, groups)
	glrt <- list(glrtReg(data, type, groups))
	names(glrt) <- "GLRT"

	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		
		GLRTpermut  <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
			p1 <- sample(x1+x2, x1, replace=FALSE)
			grps <- rep(0, x1+x2)
			grps[p1] <- 1
			glrt.mcnp <- glrtReg(data, type, grps)
			return(glrt.mcnp)
		}
		parallel::stopCluster(cl) 
	}else{
		GLRTpermut <- apply(as.matrix(1:numPerms), 1, function(x, x1, x2, data, type){	
					p1 <- sample(x1+x2, x1, replace=FALSE)
					grps <- rep(0, x1+x2)
					grps[p1] <- 1	
					glrt.mcnp <- glrtReg(data, type, grps)
					return(glrt.mcnp)
				}, x1=x1, x2=x2, data=data, type=type)
	}
	
	pvalue <- list(sum(GLRTpermut >= glrt)/numPerms)
	names(pvalue) <- "pvalue"
	
	#Label our output list
	l <- list(reg, glrt, pvalue)
	keys <- unique(unlist(lapply(l, names)))
	results	<- setNames(do.call(mapply, c(FUN=c, lapply(l, `[`, keys))), keys)	
	
	return(results)
}

graphNetworkPlot <- function(data, type, main="Network Plot", labels, groupCounts, groupLabels){
	if(missing(data) || missing(type))
		stop("data and/or type is missing.")
	
	#Set up plot dot colors
	if(missing(groupCounts)){
		myColors <- "red"
	}else{
		allColors <- rainbow(length(groupCounts)+1)
		myColors <- NULL
		for(i in 1:length(groupCounts))
			myColors <- c(myColors, rep(allColors[i], groupCounts[i]))
	}
	
	#Take only the first column of data if it is multi columned
	if(class(data) == "data.frame" || class(data) == "matrix")
		data <- data[,1]
	
	y <- vec2mat(data, type)
	g <- network::network(as.matrix(y), directed=FALSE)
	if(!missing(labels))
		network::network.vertex.names(g) <- as.data.frame(labels)
	
	network::plot.network(g, mode="circle", vertex.col=myColors, label=network::network.vertex.names(g), main=main, edge.col="black")
	
	if(!missing(groupLabels))
		legend("topright", legend=groupLabels, fill=allColors, horiz=FALSE)
}

pairedPvalue  <- function(data, type, groups, numPerms=10, parallel=FALSE, cores=3){	
	if(missing(data) || missing(type) || missing(groups))
		stop("data, type, and/or groups is missing.")
	
	if(numPerms <= 0)
		stop("The number of permutations must be an integer greater than 0.")
	
	grp1 <- data[,groups==unique(groups)[1]]
	grp2 <- data[,groups==unique(groups)[2]]
	groupSize <- ncol(grp1)
	if(groupSize != ncol(grp2))
		stop("Groups must be the same size.")
	
	data <- cbind(grp1, grp2)
	gstarDistance <- calcDistance(estGStar(grp1), estGStar(grp2), type) 
	
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		
		permDistances <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
			samps <- sample(0:1, groupSize, replace=TRUE)
			samps <- samps*groupSize + 1:groupSize
			
			gstar1 <- estGStar(data[,samps])
			gstar2 <- estGStar(data[,-samps])
			return(calcDistance(gstar1, gstar2))
		}	
		parallel::stopCluster(cl) 
	}else{
		permDistances <- NULL
		for(i in 1:numPerms){ 	
			samps <- sample(0:1, groupSize, replace=TRUE)
			samps <- samps*groupSize + 1:groupSize
			
			gstar1 <- estGStar(data[,samps])
			gstar2 <- estGStar(data[,-samps])
			
			permDistances <- c(permDistances, calcDistance(gstar1, gstar2))
		}
	}
	
	Pvalue <- sum(permDistances>=gstarDistance)/numPerms 
	return(Pvalue)
}

plotMDS <- function(data, groups, estGstar=TRUE, paired=FALSE, returnCoords=FALSE, ...){
	if(missing(data))
		stop("data is missing.")
	
	if(missing(groups))
		groups <- rep(1, ncol(data))
	
	if(length(groups) != ncol(data))
		stop("groups must be the same length as the data set.")
	
	if(min(groups) <= 0)
		groups <- groups - min(groups) + 1
	
	dlen <- ncol(data)
	uni <- unique(groups)
	mycol <- groups + 1
	mycex <- rep(1, length(groups))
	numcols <- max(length(uni), 2, max(mycol))
	
	if(paired){
		if(length(uni) != 2 || sum(groups==uni[1]) != sum(groups==uni[2]))
			stop("There must be exactly TWO groups of the same length when using paired==TRUE.")
		groupSize <- ncol(data)/2
	}
	
	if(estGstar){
		gss <- NULL
		if(length(uni) > 1){
			for(i in 1:length(uni)){
				if(sum(groups==uni[i]) == 1)
					next
				gss <- cbind(gss, estGStar(data[,groups==uni[i]]))
				colnames(gss)[ncol(gss)] <- paste("Group", i, "Gstar")
				mycol <- c(mycol, max(mycol)+1)
			}
		}
		gss <- cbind(gss, estGStar(data))
		numcols <- numcols + ncol(gss)
		colnames(gss)[ncol(gss)] <- paste("Combined Gstar")
		mycol <- c(mycol, 1)
		mycex <- c(mycex, rep(1.5, ncol(gss)))
		data <- cbind(data, gss)
	}
	palette(c("black", rainbow(numcols)))
	
	loc <- cmdscale(dist(t(data), method="binary"), k=2) 
	plot(loc, xlab="MDS 1", ylab="MDS 2", pch=16, col=mycol, cex=mycex, ...)
	
	if(paired){
		for(i in 1:groupSize)						
			segments(loc[i,1], loc[i,2], loc[i+groupSize,1], loc[i+groupSize,2], lty=3)
		
		if(estGstar)
			segments(loc[groupSize*2+1,1], loc[groupSize*2+1,2], loc[groupSize*2+2,1], loc[groupSize*2+2,2], lty=3) 	
	}
	
	if(estGstar){
		ptr <- nrow(loc) - ncol(gss)
		for(i in 1:ncol(gss)){
			text(loc[ptr+i,1], loc[ptr+i,2], rownames(loc)[ptr+i], pos=3)
		}
	}
	
	if(returnCoords)
		return(loc)
}

plotHeatmap <- function(data, type, names, ...){
	if(missing(data) || missing(type))
		stop("data and/or type is missing.")
	
	#Take only the first column of data if it is multi columned
	if(class(data) == "data.frame" || class(data) == "matrix")
		data <- data[,1]
	
	mat <- vec2mat(data, type)
	
	if(missing(names))
		names <- 1:ncol(mat)
	colnames(mat) <- names
	rownames(mat) <- names
	
	colfunc <- colorRampPalette(c("blue", "grey"))
	
	gplots::heatmap.2(mat, symm=TRUE, Rowv=NA, dendrogram="none", trace="none", col=colfunc(10), ...)
}

lrtPvalue <- function(data, type, groups, numPerms=10, parallel=FALSE, cores=3){
	if(missing(data) || missing(type) || missing(groups))
		stop("data, type, and/or groups is missing.")
	
	if(numPerms <= 0)
		stop("The number of permutations must be an integer greater than 0.")
	
	if(length(unique(groups)) != 2)
		stop("There must be exactly two groups.")
	
	x1 <- sum(groups==unique(groups)[1])
	x2 <- sum(groups==unique(groups)[2])
	if(x1 <= 0 || x2 <= 0)
		stop("Each group must have at least 1 subject.")
	
	grp1 <- data[,groups==unique(groups)[1]]
	grp2 <- data[,groups==unique(groups)[2]]
	groupSize <- ncol(grp1)
	
	gstar <- estGStar(data)
	gstar1 <- estGStar(grp1)
	gstar2 <- estGStar(grp2)
	
	tau <- estTau(data, type, gstar)
	tau1 <- estTau(grp1, type, gstar1)
	tau2 <- estTau(grp2, type, gstar2)
	
	ll <- estLogLik(data, type, gstar, tau)
	ll1 <- estLogLik(grp1, type, gstar1, tau1)
	ll2 <- estLogLik(grp2, type, gstar2, tau2)
	
	e10raw <- -2 * (ll-ll1-ll2)
	
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		
		lambda  <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
			g <- sample(1:ncol(data), groupSize)
			grp1 <- data[,g]
			grp2 <- data[,-g]
			
			gstar1 <- estGStar(grp1)
			gstar2 <- estGStar(grp2)
			tau1 <- estTau(grp1, type, gstar1)
			tau2 <- estTau(grp2, type, gstar2)
			ll1 <- estLogLik(grp1, type, gstar1, tau1)
			ll2 <- estLogLik(grp2, type, gstar2, tau2)
			
			e10 <- -2 * (ll-ll1-ll2)
			return(e10)
		}
		parallel::stopCluster(cl) 
	}else{
		lambda <- NULL
		for(i in 1:numPerms){
			g <- sample(1:ncol(data), groupSize)
			grp1 <- data[,g]
			grp2 <- data[,-g]
			
			gstar1 <- estGStar(grp1)
			gstar2 <- estGStar(grp2)
			tau1 <- estTau(grp1, type, gstar1)
			tau2 <- estTau(grp2, type, gstar2)
			ll1 <- estLogLik(grp1, type, gstar1, tau1)
			ll2 <- estLogLik(grp2, type, gstar2, tau2)
			
			e10 <- -2 * (ll-ll1-ll2)
			lambda <- c(lambda, e10)
		}
	}
	pValue <- sum(lambda > e10raw)/numPerms	
	
	return(pValue)
}

gaConsensus <- function(data, groups, iters=10, nRuns=1, popSize=200, method="manhattan", parallel=FALSE, cores=3){
	if(missing(data) || missing(groups))
		stop("data and/or groups is missing.")
	
	if(length(unique(groups)) != 2)
		stop("There must be exactly two groups.")
	
	### Get our groups in the right order
	x1 <- sum(groups==unique(groups)[1])
	x2 <- sum(groups==unique(groups)[2])
	if(x1 <= 0 || x2 <= 0)
		stop("Each group must have at least 1 subject.")
	
	grp1 <- data[,groups==unique(groups)[1]]
	grp2 <- data[,groups==unique(groups)[2]]
	groupSize1 <- ncol(grp1)
	groupSize2 <- ncol(grp2)
	data <- cbind(grp1, grp2)
	
	### Set up our function for the GA
	b <- dist(c(rep(0, groupSize1), rep(1, groupSize2)), method)
	myScore <- function(indices) {
		if(sum(indices) == 0)
			return(0)
		edges <- which(indices==1)
		
		a <- Reduce("+", dists[edges])
		mycor <- cor(a, b)
		return(-mycor)
	}
	
	### Precompute all our distances
	dists <- vector("list", nrow(data))
	for(i in 1:nrow(data))
		dists[[i]] <- dist(as.vector(t(data[i,])), method)
	
	### Run the GA
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		
		res <- foreach::foreach(i=1:nRuns, .combine=rbind, .inorder=FALSE, .multicombine=TRUE, .packages=c("genalg")) %dopar%{
			tempGa <- genalg::rbga.bin(size=nrow(data), iters=iters, popSize=popSize, evalFunc=myScore, verbose=TRUE)
			return(tempGa$population[which(tempGa$evaluation==min(tempGa$evaluation)),])
		}
		parallel::stopCluster(cl) 
	}else{
		res <- NULL
		for(i in 1:nRuns){
			tempGa <- genalg::rbga.bin(size=nrow(data), iters=iters, popSize=popSize, evalFunc=myScore, verbose=TRUE)	
			res <- rbind(res, tempGa$population[which(tempGa$evaluation==min(tempGa$evaluation)),])
		}
	}
	
	### Clean up our results
	res <- t(res)
	res <- res[,!duplicated(t(res)), drop=FALSE]
	corrs <- apply(res, 2, myScore)
	res <- res[,order(corrs), drop=FALSE]
	
	colnames(res) <- paste("Solution", 1:ncol(res))
	corrs <- corrs[order(corrs)] * -1
	
	return(list(solutions=res, corrs=corrs))
}


