#' remakeRandomData
#' 
#' Given an object output by \code{makeRandomData}, make another data set with the same
#' distribution. 
#' 
#' @param n A \code{numeric} sample size
#' @param object An object of class \code{"makeRandomData"}
#' @param setA Value(s) to set treatment variable to. If \code{NULL} then treatment
#' is simulated according to observed data functions. If \code{length(A)==1} then all
#' values are set to this single value. Otherwise, if A is a vector of length n, A is set
#' according to this vector. 
#' @param setW A matrix of proper size with values of covariates set to fixed levels. Useful for
#' plotting methods. Assumes proper dimensions. If \code{NULL} simulates W according to distributions.
#' @return An object of class \code{"makeRandomData"} with the following entries
#' \item{W}{A matrix of covariates}
#' \item{A}{A vector of binary treatments}
#' \item{Y}{A vector of continuously valued outcome}
#' \item{distW}{A list containing relevant information needed to reproduce data sets}
#' \item{fnG0}{A list of lists containing relevant information needed to reproduce data sets}
#' \item{fnQ0}{A list of lists containing relevant information needed to reproduce data sets}
#' \item{distErrY}{A list containing relevant information needed to reproduce data sets}
#' @export

remakeRandomData <- function(n, object, setA = NULL, setW = NULL, ...){

	# draw random number of covariates
	D <- ncol(object$W)

	#----------------------------------------------------------------------
	# Simulate W	
	#----------------------------------------------------------------------
	# initialize empties	
	W <- matrix(nrow = n, ncol = D)
	distW <- vector(mode = "list", length = D)
	if(is.null(setW)){
		for(d in 1:D){
			if(d == 1){
				W[,d] <- do.call(object$distW[[d]]$fn, args = c(list(n=n), object$distW[[d]]$parm))			
			}else{
				W[,d] <- do.call(object$distW[[d]]$fn, args = c(list(n=n, x = W[,d-1]), object$distW[[d]]$parm))
			}
		}
	}else{
		W <- as.matrix(setW)
	}

	#----------------------------------------------------------------------
	# Simulate propensity -- only if setA == NULL
	#----------------------------------------------------------------------
	if(is.null(setA)){
		# draw random number of main terms
		Mg1 <- length(object$fnG0$uni)
		Mg2 <- length(object$fnG0$biv)
		Mg3 <- length(object$fnG0$tri)
		Mg4 <- length(object$fnG0$quad)

		# initialize empty
		logitg0 <- rep(0, n)

		# univariate
		if(Mg1 > 0){
			for(m in 1:Mg1){
				fOut <- do.call(object$fnG0$uni[[m]]$fn, args = c(list(x = W[,object$fnG0$uni[[m]]$whichColsW]), object$fnG0$uni[[m]]$parm))
				# add to current logitg0
				logitg0 <- logitg0 + fOut
			}
		}

		# two-way interactions
		if(Mg2 > 0){
			for(m in 1:Mg2){
				# call function with parameters
				fOut <- do.call(object$fnG0$biv[[m]]$fn, 
				                args = c(list(x1 = W[,object$fnG0$biv[[m]]$whichColsW[1]],
				                              x2 = W[,object$fnG0$biv[[m]]$whichColsW[2]]),
				                object$fnG0$biv[[m]]$parm))
				# add to current logitg0
				logitg0 <- logitg0 + fOut
			}
		}

		#trivariate
		if(Mg3 > 0){
			for(m in 1:Mg3){
				# call function with parameters
				fOut <- do.call(object$fnG0$tri[[m]]$fn, 
				                args = c(list(x1 = W[,object$fnG0$tri[[m]]$whichColsW[1]],
				                              x2 = W[,object$fnG0$tri[[m]]$whichColsW[2]],
				                              x3 = W[,object$fnG0$tri[[m]]$whichColsW[3]]),
				                object$fnG0$tri[[m]]$parm))# save output in list
				# add to current logitg0
				logitg0 <- logitg0 + fOut
			}
		}	
		#quadravariant	
		if(Mg4 > 0){
			for(m in 1:Mg4){
				# call function with parameters
				fOut <- do.call(object$fnG0$quad[[m]]$fn, 
				                args = c(list(x1 = W[,object$fnG0$quad[[m]]$whichColsW[1]],
				                              x2 = W[,object$fnG0$quad[[m]]$whichColsW[2]],
				                              x3 = W[,object$fnG0$quad[[m]]$whichColsW[3]],
				                              x4 = W[,object$fnG0$quad[[m]]$whichColsW[4]]),
				                object$fnG0$quad[[m]]$parm))# save output in list
				# add to current logitg0
				logitg0 <- logitg0 + fOut
			}
		}

   	  logitg0 <- .8^object$its*logitg0 + .8^object$its*object$skewage

		# correct for positivity violations
		logitg0[plogis(logitg0) < object$minG0] <- qlogis(object$minG0)
		logitg0[plogis(logitg0) > 1 - object$minG0] <- qlogis(1 - object$minG0)
		
		# simulate A
		A <- rbinom(n, 1, plogis(logitg0))
	}else{
		if(length(setA)==1){
			A <- rep(setA, n)			
		}else{
			A <- setA
		}
		logitg0 <- Inf
	}
	# matrix with A and W
	AW <- cbind(A, W)

	#----------------------------------------------------------------------
	# Simulate Y
	#----------------------------------------------------------------------
	# draw random number of main terms between 2 and D + 1, where we set 
	# 2 to be the minimum so that we ensure there is confounding. 
	MQ1 <- length(object$fnQ0$uni)
	MQ2 <- length(object$fnQ0$biv)
	MQ3 <- length(object$fnQ0$tri)
	MQ4 <- length(object$fnQ0$quad)
	
	# empty 
	Q0 <- rep(0, n)

	# main terms
	if(MQ1 > 0){
		for(m in 1:MQ1){
			fOut <- do.call(object$fnQ0$uni[[m]]$fn, args = c(list(x = AW[,object$fnQ0$uni[[m]]$whichColsAW]), object$fnQ0$uni[[m]]$parm))
			# add 
			Q0 <- Q0 + fOut
		}
	}
	# two-way interactions
	if(MQ2 > 0){
		for(m in 1:MQ2){
			# call function with parameters 
			fOut <- do.call(object$fnQ0$biv[[m]]$fn, 
			                args = c(list(x1 = AW[,object$fnQ0$biv[[m]]$whichColsAW[1]], 
			                              x2 = AW[,object$fnQ0$biv[[m]]$whichColsAW[2]]),
			                object$fnQ0$biv[[m]]$parm))
			# add
			Q0 <- Q0 + fOut
		}
	}

	# three-way interactions 
	if(MQ3 > 0){
		for(m in 1:MQ3){
			# call function with parameters
			fOut <- do.call(object$fnQ0$tri[[m]]$fn, 
			                args = c(list(x1 = AW[,object$fnQ0$tri[[m]]$whichColsAW[1]], 
			                              x2 = AW[,object$fnQ0$tri[[m]]$whichColsAW[2]],
			                              x3 = AW[,object$fnQ0$tri[[m]]$whichColsAW[3]]),
			                object$fnQ0$tri[[m]]$parm))			
			# add
			Q0 <- Q0 + fOut
		}
	}
	# four-way interactions 
	if(MQ4 > 0){
		for(m in 1:MQ4){
			# call function with parameters
			fOut <- do.call(object$fnQ0$quad[[m]]$fn, 
			                args = c(list(x1 = AW[,object$fnQ0$quad[[m]]$whichColsAW[1]], 
			                              x2 = AW[,object$fnQ0$quad[[m]]$whichColsAW[2]],
			                              x3 = AW[,object$fnQ0$quad[[m]]$whichColsAW[3]],
			                              x4 = AW[,object$fnQ0$quad[[m]]$whichColsAW[4]]),
			                object$fnQ0$quad[[m]]$parm))			
			# add
			Q0 <- Q0 + fOut
		}
	}


	# Drawing an error function
	# evaluate error function
	errOut <- do.call(object$distErrY$fn, args = c(list(AW=AW, n=n), object$distErrY$parm))

	# compute Y
	Y <- Q0 + errOut * object$errMult

	return(list(W = W, A = A, Y = Y, Q0 = Q0, g0 = plogis(logitg0), err = errOut))
}