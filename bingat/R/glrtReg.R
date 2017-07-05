glrtReg <-
function(data, type, groups){
	regCoeffHA <- glmReg(data, type, groups)
	regCoeffH0 <- glmReg(data, type, 0)
	
	glrt <- 2 * (regCoeffHA$loglik - regCoeffH0$loglik)
	return(glrt)
}
