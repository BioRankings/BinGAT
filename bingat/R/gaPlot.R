gaPlot <-
function(evalSumm){
	graphics::plot(evalSumm[,4], type="l", ylab="Score", ylim=c(0, 1), lwd=2, main="Eval Scores by Iteration", xlab="Iteration")
	graphics::lines(evalSumm[,6], col="red", lwd=2)
	graphics::lines(evalSumm[,1], col="blue", lwd=2)
	graphics::legend("topleft", colnames(evalSumm)[c(4, 6, 1)], pch=16, col=c("black", "red", "blue"))
}
