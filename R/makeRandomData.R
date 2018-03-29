#' makeRandomData
#'
#' Simulate a random data generating distribution. See Details to see how it's done.
#'
#'
#' Draw random data generating distribution. The code simulates as follows: \cr
#' Simulate W \itemize{
#' \item Choose \code{D}, a random number of covariates from Uniform(1,maxD)
#' \item Randomly choose \code{d1}, the distribution for first covariate, from \code{distW.uncor}
#' \item Call \code{paste0(d1,".Parm")} (which is assumed to exist in the global environment) to get parameters for \code{d1}
#' \item Call \code{d1} with parameters chosen in previous step to generate first covariate
#' \item Randomly choose \code{d2}, the distribution for second covariate, from \code{c(distW.uncor, distW.cor)} (now including correlated distributions).
#' \item Call \code{paste0(d2,".Parm")} to get parameters for \code{d2}
#' \item Call \code{d2} with parameters chosen in previous step and first covariate (to possibly induce correlation with previous covariate) to generate second covariate.
#' \item Continue until \code{maxD} is reached.
#' }
#'
#' Simulating A \itemize{
#' \item Draw \code{MG1}, a random number of main terms for propensity score, from Uniform(1,maxD)
#' \item For \code{i} = 1,...,\code{MG1}, draw \code{f} a random function from \code{funcG0.uni}
#' \item Call \code{paste0(f,"Parm")} (assumed to exist) to generate parameters for \code{f}
#' \item Call \code{f} with \code{x = W[,i]} to generate \code{i}-th main term.
#' \item If \code{D > 1}, draw \code{MG2}, a random number of bivariate interactions for propensity, from Uniform(1,\code{MG1-1})
#' \item For \code{i} = 1,...,\code{MG2}, draw \code{f} a random function from \code{funcG0.biv}
#' \item Call \code{paste0(f,"Parm")} (assumed to exist) to generate parameters for \code{f}
#' \item Draw \code{j} and \code{k}, two random columns of \code{W}
#' \item Call \code{f} with \code{x1 = W[,j], x2 = W[,k]} to generate \code{i}-th bivariate interaction term.
#' \item If \code{D > 2}, draw \code{MG3}, a random number of trivariate interactions for propensity, from Uniform(1,\code{MG2-1})
#' \item For \code{i} = 1,...,\code{MG3}, draw \code{f} a random function from \code{funcG0.tri}
#' \item Call \code{paste0(f,"Parm")} (assumed to exist) to generate parameters for \code{f}
#' \item Draw \code{j},\code{k},\code{l}, three random columns of \code{W}
#' \item Call \code{f} with \code{x1 = W[,j], x2 = W[,k], x3 = W[,l]} to generate \code{i}-th trivariate interaction term.
#' \item Sum together all main, bivariate, and trivariate terms to get the logit propensity score.
#' \item Divide the logit propensity score by 1.01 until the minimum propensity score is bigger than \code{minG0}.
#' \item Draw A from Bernoulli distribution with conditional probability equal to the adjusted propensity score.
#' }
#'
#' Simulating Y \itemize{
#' \item Draw \code{MQ1}, a random number of main terms for outcome regression score, from Uniform(2,maxD)
#' \item For \code{i} = 1,...,\code{MQ1}, draw \code{f} a random function from \code{funcQ0.uni}
#' \item Call \code{paste0(f,"Parm")} (assumed to exist) to generate parameters for \code{f}
#' \item Call \code{f} with \code{x = AW[,i]} to generate \code{i}-th main term, where AW is matrix with first column A, other columns W.
#' \item If \code{D > 1}, draw \code{MQ2}, a random number of bivariate interactions for outcome regression, from Uniform(1,\code{MQ1-1})
#' \item For \code{i} = 1,...,\code{MQ2}, draw \code{f} a random function from \code{funcQ0.biv}
#' \item Call \code{paste0(f,"Parm")} (assumed to exist) to generate parameters for \code{f}
#' \item Draw \code{j} and \code{k}, two random columns of \code{W}
#' \item Call \code{f} with \code{x1 = W[,j], x2 = W[,k]} to generate \code{i}-th bivariate interaction term.
#' \item If \code{D > 2}, draw \code{MQ3}, a random number of trivariate interactions for propensity, from Uniform(1,\code{MQ2-1})
#' \item For \code{i} = 1,...,\code{MQ3}, draw \code{f} a random function from \code{funcQ0.tri}
#' \item Call \code{paste0(f,"Parm")} (assumed to exist) to generate parameters for \code{f}
#' \item Draw \code{j},\code{k},\code{l}, three random columns of \code{W}
#' \item Call \code{f} with \code{x1 = W[,j], x2 = W[,k], x3 = W[,l]} to generate \code{i}-th trivariate interaction term.
#' \item Sum together all main, bivariate, and trivariate terms to get the conditional mean outcome, \code{Q0}.
#' \item Draw \code{errf}, arandom function from \code{errY} from which to draw errors
#' \item Call \code{paste0(errf,"Parm"} (assumed to exist) to generate random parameters for error distribution.
#' \item Call \code{errf} with randomly chosen parameters to draw \code{e}, an error for the outcome.
#' \item Let \code{Y = Q0 + e}.
#' }
#'
#' @param n A \code{numeric} sample size
#' @param maxD A \code{numeric} indicating maximum number of covariates
#' @param distW.uncor A \code{vector} of \code{characters} that are functions in the global environment that
#' generate covariates. These functions should have an associated parameter function with \code{"Parm"} added
#' to the end of their name. See the format of the functions provided to understand their structure.
#' @param distW.cor Same as above, but generates a correlated covariate
#' @param funcG0.uni A \code{vector} of \code{characters} that are functions in the global environment that
#' are used to generate main terms for logit of the propensity score.
#' @param funcG0.biv Ditto above, but for bivariate functions.
#' @param funcG0.tri Ditto above, but for trivariate functions.
#' @param funcQ0.uni Ditto \code{funcG0.uni}, but for the outcome regression
#' @param funcQ0.biv Ditto
#' @param funcQ0.tri Ditto
#' @param errY A \code{vector} of \code{characters} that are functions in the global environment that
#' are used to generate error terms for the outcome regression.
#' @param minG0 The minimum value for the propensity score (default 0.01).
#'
#' @return An object of class \code{"makeRandomData"} with the following entries
#' \item{W}{A matrix of covariates}
#' \item{A}{A vector of binary treatments}
#' \item{Y}{A vector of continuously valued outcome}
#' \item{distW}{A list containing relevant information needed to reproduce data sets}
#' \item{fnG0}{A list of lists containing relevant information needed to reproduce data sets}
#' \item{fnQ0}{A list of lists containing relevant information needed to reproduce data sets}
#' \item{distErrY}{A list containing relevant information needed to reproduce data sets}
#' \item{divideLogitG0}{A numeric of the scaling factor for the propensity}
#' @export

makeRandomData <- function(n,
                           maxD,
                           minObsA = 30,
                           func.distW = c("uniformW","normalW","bernoulliW","binomialW","gammaW",
                                     "normalWCor","bernoulliWCor","uniformWCor", "gammaPointMassW",
                                     "binomialFracW","normalPointMassW"),
                           funcG0.uni = c("linUni","polyUni","sinUni","jumpUni","pLogisUni",
                                          "dNormUni","qGammaUni","dNormMixUni","cubicSplineUni",
                                          "linSplineUni"),
                           funcG0.biv = c("linBiv","polyBiv","sinBiv","jumpBiv",
                               "dNormAddBiv","dNormMultBiv","polyJumpBiv",
                               "sinJumpBiv","dNormMultJumpBiv","linSplineBiv"),
                           funcG0.tri = c("linTri", "polyTri", "sinTri", "jumpTri",
	                            "linJump2Tri","linJump1Tri","polyJump2Tri",
	                            "polyJump1Tri","sinJump2Tri","sinJump1Tri",
	                            "dNormAddTri","dNormMultTri","dNormMultJump2Tri",
	                            "dNormMultJump1Tri","pLogisAddTri"),
                           funcG0.quad = c("linQuad", "polyQuad", "sinQuad", "jumpQuad",
	                            "linJump2Quad","linJump1Quad","polyJump2Quad",
	                            "polyJump1Quad","sinJump2Quad","sinJump1Quad",
	                            "dNormAddQuad","dNormMultQuad","dNormMultJump2Quad",
	                            "dNormMultJump1Quad","pLogisAddQuad"),
                           funcQ0.uni = c("linUni","polyUni","sinUni","jumpUni","pLogisUni",
                                          "dNormUni","qGammaUni","dNormMixUni","cubicSplineUni",
                                          "linSplineUni"),
                           funcQ0.biv = c("linBiv","polyBiv","sinBiv","jumpBiv",
                               "dNormAddBiv","dNormMultBiv","polyJumpBiv",
                               "sinJumpBiv","dNormMultJumpBiv","linSplineBiv"),
                           funcQ0.tri = c("linTri", "polyTri", "sinTri", "jumpTri",
	                            "linJump2Tri","linJump1Tri","polyJump2Tri",
	                            "polyJump1Tri","sinJump2Tri","sinJump1Tri",
	                            "dNormAddTri","dNormMultTri","dNormMultJump2Tri",
	                            "dNormMultJump1Tri","pLogisAddTri"),
                           funcQ0.quad = c("linQuad", "polyQuad", "sinQuad", "jumpQuad",
	                            "linJump2Quad","linJump1Quad","polyJump2Quad",
	                            "polyJump1Quad","sinJump2Quad","sinJump1Quad",
	                            "dNormAddQuad","dNormMultQuad","dNormMultJump2Quad",
	                            "dNormMultJump1Quad","pLogisAddQuad"),
                           errY = c("normalErr","uniformErr","gammaErr","normalErrAW","uniformErrAW"),
                           minG0 = 1e-3,
                           minR2 = 0.01, maxR2 = 0.99,
                           ...){
	# draw random number of covariates
  D <- sample(1:maxD, 1)
	#----------------------------------------------------------------------
	# Simulate W
	#----------------------------------------------------------------------
	# initialize empties
	W <- matrix(nrow = n, ncol = D)
	distW <- vector(mode = "list", length = D)
	for(d in 1:D){
		# randomly sample distribution
		if(d==1){
			# sample distribution
			thisD <- sample(func.distW, 1)
			# generate parameters
			thisParm <- do.call(paste0(thisD,"Parm"), args = list(W = NULL))
			# generate covariate
			W[,d] <- do.call(thisD, args = c(list(n=n, W = NULL), thisParm))
		}else{
			# sample distribution
			thisD <- sample(func.distW, 1)
			# generate parameters
			thisParm <- do.call(paste0(thisD,"Parm"), args = list(W = data.frame(W[,1:(d-1),drop=FALSE])))
			# generate covariate
			W[,d] <- do.call(thisD, args = c(list(n=n, W = data.frame(W[,1:(d-1),drop=FALSE])), thisParm))
		}
		# save distributions
		distW[[d]]$fn <- thisD; distW[[d]]$parm <- thisParm
	}


	#----------------------------------------------------------------------
	# Simulate propensity
	#----------------------------------------------------------------------
	# draw random number of main terms
	Mg1 <- round(runif(1, -0.5, D + 0.5))

	# draw random number of two-way interaction terms
	Mg2 <- round(runif(1, -0.5, D + 0.5))

	# draw random number of three-way interaction terms
	Mg3 <- round(runif(1, -0.5, D + 0.5))

	# draw random number of four-way interaction terms
	Mg4 <- round(runif(1, -0.5, D + 0.5))

	# initialize empty
	logitg0 <- rep(0, n)

	# make sure there are at least some A in each group
	A <- rep(0, n)
	while(sum(A == 1) < minObsA | sum(A==0) < minObsA){
		# univariate
		if(Mg1 > 0){
			uniG0 <- vector(mode="list", length = Mg1)
			for(m in 1:Mg1){
				# draw random function
				thisF <- sample(funcG0.uni, 1)
				# draw random column of W
				wCol <- sample(1:ncol(W), 1)
				# get parameters
				thisParm <- do.call(paste0(thisF,"Parm"), args = list(x = W[,wCol]))
				# call function with parameters
				fOut <- do.call(thisF, args = c(list(x = W[,wCol]), thisParm))
				# save output in list
				uniG0[[m]] <- list(fn = thisF, parm = thisParm, whichColsW = wCol)
				# add to current logitg0
				logitg0 <- logitg0 + fOut
			}
		}else{
			uniG0 <- NULL
		}
		# two-way interactions
		if(Mg2 > 0 &  D > 1){
			# empty list
			bivG0 <- vector(mode="list",length=Mg2)
			# all two-way column combinations
			comb <- as.matrix(combn(D,2))
			# randomly sample Mg2 two-way interactions without replacement
			combCols <- sample(1:ncol(comb),Mg2,replace = TRUE)
			for(m in 1:Mg2){
				# the two columns used for this interaction
				theseCols <- comb[,combCols[m]]
				# the random function to be used
				thisF <- sample(funcG0.biv, 1)
				# get parameters for function
				thisParm <- do.call(paste0(thisF,"Parm"), args = list(x1 = W[,theseCols[1]], x2 = W[,theseCols[2]]))
				# call function with parameters
				fOut <- do.call(thisF, args = c(list(x1 = W[,theseCols[1]], x2 = W[,theseCols[2]]),thisParm))
				# save output in list
				bivG0[[m]] <- list(fn = thisF, parm = thisParm, whichColsW = theseCols)
				# add to current logitg0
				logitg0 <- logitg0 + fOut
			}
		}else{
			bivG0 <- NULL
		}

		#trivariate
		if(Mg3 > 0 & D > 2){
			# empty list
			triG0 <- vector(mode="list",length=Mg3)
			# all three way choices of columns
			comb <- as.matrix(combn(D, 3))
			# randomly sample Mg3 three-way interactions without replacement
			combCols <- sample(1:ncol(comb),Mg3,replace = TRUE)
			for(m in 1:Mg3){
				# the three columns used for this function
				theseCols <- comb[,combCols[m]]
				# the random function to be used
				thisF <- sample(funcG0.tri, 1)
				# get parameters for function
				thisParm <- do.call(paste0(thisF,"Parm"), args = list(x1 = W[,theseCols[1]], x2 = W[,theseCols[2]], x3=W[,theseCols[3]]))
				# call function with parameters
				fOut <- do.call(thisF, args = c(list(x1 = W[,theseCols[1]], x2 = W[,theseCols[2]], x3=W[,theseCols[3]]),thisParm))
				# save output in list
				triG0[[m]] <- list(fn = thisF, parm = thisParm, whichColsW = theseCols)
				# add to current logitg0
				logitg0 <- logitg0 + fOut
			}
		}else{
			triG0 <- NULL
		}
		# quadravariate
		if(Mg4 > 0 & D > 3){
			# empty list
			quadG0 <- vector(mode="list",length=Mg4)
			# all four way choices of columns
			comb <- as.matrix(combn(D, 4))
			# randomly sample Mg3 four-way interactions without replacement
			combCols <- sample(1:ncol(comb), Mg4, replace = TRUE)
			for(m in 1:Mg4){
				# the four columns used for this function
				theseCols <- comb[,combCols[m]]
				# the random function to be used
				thisF <- sample(funcG0.quad, 1)
				# get parameters for function
				thisParm <- do.call(paste0(thisF,"Parm"),
				                    args = list(x1 = W[,theseCols[1]],
				                    x2 = W[,theseCols[2]], x3=W[,theseCols[3]],
				                    x4 = W[,theseCols[4]]))
				# call function with parameters
				fOut <- do.call(thisF, args = c(list(x1 = W[,theseCols[1]],
				                                     x2 = W[,theseCols[2]],
				                                     x3 = W[,theseCols[3]],
				                                     x4 = W[,theseCols[4]]),thisParm))
				# save output in list
				quadG0[[m]] <- list(fn = thisF, parm = thisParm, whichColsW = theseCols)
				# add to current logitg0
				logitg0 <- logitg0 + fOut
			}
		}else{
			quadG0 <- NULL
		}

		# correct for positivity violations
		logitg0[plogis(logitg0) < minG0] <- qlogis(minG0)
		logitg0[plogis(logitg0) > 1 - minG0] <- qlogis(1 - minG0)

		# simulate A
		A <- rbinom(n, 1, plogis(logitg0))
	}
	# matrix with A and W
	AW <- cbind(A, W)

	#----------------------------------------------------------------------
	# Simulate Y
	#----------------------------------------------------------------------
	# draw random number of main terms between 0 and D + 1, where we set
	# no minimum, which allows there to be cases of just pure noise
	MQ1 <- round(runif(1, -0.5, D + 1.5))

	# draw random number of interaction terms
	MQ2 <- round(runif(1, -0.5, D + 1.5))
	MQ3 <- round(runif(1, -0.5, D + 1.5))
	MQ4 <- round(runif(1, -0.5, D + 1.5))

	# empty
	Q0 <- rep(0, n)

	# main terms
	# empty
	if(MQ1 > 0){
		uniQ0 <- vector(mode="list", length = MQ1)
		for(m in 1:MQ1){
			# randomly draw function
			thisF <- sample(funcQ0.uni, 1)
			# randomly draw column of AW
			awCol <- sample(1:ncol(AW), 1)
			# get parameters
			thisParm <- do.call(paste0(thisF,"Parm"), args = list(x = AW[,awCol]))
			# call function with parameters
			fOut <- do.call(thisF, args = c(list(x = AW[,awCol]), thisParm))
			# save results
			uniQ0[[m]] <- list(fn = thisF, parm = thisParm, whichColsAW = awCol)
			# add
			Q0 <- Q0 + fOut
		}
	}else{
		uniQ0 <- NULL
	}
	# two-way interactions
	if(MQ2 > 0 & D > 0){
		bivQ0 <- vector(mode="list",length=MQ2)
		# all combinations of columns
		comb <- as.matrix(combn(D+1,2))
		# randomly sample columns
		combCols <- sample(1:ncol(comb),MQ2, replace = TRUE)
		for(m in 1:MQ2){
			# what columns to use for this interaction
			theseCols <- comb[,combCols[m]]
			# randomly sample function
			thisF <- sample(funcQ0.biv,1)
			# get parameters
			thisParm <- do.call(paste0(thisF,"Parm"), args = list(x1 = AW[,theseCols[1]], x2 = AW[,theseCols[2]]))
			# call function with parameters
			fOut <- do.call(thisF, args = c(list(x1 = AW[,theseCols[1]], x2 = AW[,theseCols[2]]),thisParm))
			# save output
			bivQ0[[m]] <- list(fn = thisF, parm = thisParm, whichColsAW = theseCols)
			# add
			Q0 <- Q0 + fOut
		}
	}else{
		bivQ0 <- NULL
	}

	# three-way interactions
	if(MQ3 > 0 & D > 1){
		# empty
		triQ0 <- vector(mode="list",length=MQ3)
		# all three-way column combinations
		comb <- as.matrix(combn(D+1,3))
		# randomly sample three choices of combinations without replacement
		combCols <- sample(1:ncol(comb),MQ3, replace = TRUE)
		for(m in 1:MQ3){
			# columns to use for this interaction
			theseCols <- comb[,combCols[m]]
			# randomly sample function
			thisF <- sample(funcQ0.tri,1)
			# get parameters
			thisParm <- do.call(paste0(thisF,"Parm"), args = list(x1 = AW[,theseCols[1]], x2 = AW[,theseCols[2]], x3=AW[,theseCols[3]]))
			# call function with parameters
			fOut <- do.call(thisF, args = c(list(x1 = AW[,theseCols[1]], x2 = AW[,theseCols[2]], x3=AW[,theseCols[3]]),thisParm))
			# save output
			triQ0[[m]] <- list(fn = thisF, parm = thisParm, whichColsAW = theseCols)
			# add
			Q0 <- Q0 + fOut
		}
	}else{
		triQ0 <- NULL
	}
	# four-way interactions
	if(MQ4 > 0 & D > 2){
		# empty
		quadQ0 <- vector(mode="list",length=MQ4)
		# all four-way column combinations
		comb <- as.matrix(combn(D+1, 4))
		# randomly sample four choices of combinations without replacement
		combCols <- sample(1:ncol(comb), MQ4, replace = TRUE)
		for(m in 1:MQ4){
			# columns to use for this interaction
			theseCols <- comb[,combCols[m]]
			# randomly sample function
			thisF <- sample(funcQ0.quad,1)
			# get parameters
			thisParm <- do.call(paste0(thisF,"Parm"),
			                    args = list(x1 = AW[,theseCols[1]],
			                                x2 = AW[,theseCols[2]],
			                                x3=AW[,theseCols[3]],
			                                x4=AW[,theseCols[4]]))
			# call function with parameters
			fOut <- do.call(thisF, args = c(list(x1 = AW[,theseCols[1]],
			                                     x2 = AW[,theseCols[2]],
			                                     x3=AW[,theseCols[3]],
			                                     x4=AW[,theseCols[4]]),thisParm))
			# save output
			quadQ0[[m]] <- list(fn = thisF, parm = thisParm, whichColsAW = theseCols)
			# add
			Q0 <- Q0 + fOut
		}
	}else{
		quadQ0 <- NULL
	}


	# Drawing an error function
	# empty results
	errYList <- vector(mode = "list")
	# randomly draw error function
	errFn <- sample(errY, 1)
	# get parameters of error function
	errParm <- do.call(paste0(errFn, "Parm"), args = list(AW = AW))
	# evaluate error function
	errOut <- do.call(errFn, args = c(list(AW=AW, n=n), errParm))
	# save output
	errYList$fn <- errFn; errYList$parm <- errParm

	# make sure R2 is falling in proper range
	ct <- 0
	currR2 <- Inf
	while(currR2 < minR2 | currR2 > maxR2){
		ct <- ct + 1
		if(currR2 > maxR2){
			mult <- 1.1^ct
		}else if(currR2 < minR2){
			mult <- 1/1.1^ct
		}
		# compute Y
		Y <- Q0 + errOut
		out <- list(W = W, A = data.frame(A=A), Y = data.frame(Y=Y), distW = distW,
					minG0 = minG0, minR2 = minR2, maxR2 = maxR2,
		            fnG0 = list(uni = uniG0, biv = bivG0, tri = triG0, quad = quadG0),
			        fnQ0 = list(uni = uniQ0, biv = bivQ0, tri = triQ0, quad = quadQ0), distErrY = errYList,
			        Q0 = Q0, g0 = plogis(logitg0), errMult = mult, funcG0.uni = funcG0.uni, funcG0.biv = funcG0.biv, funcG0.tri = funcG0.tri, funcG0.quad = funcG0.quad,
			        funcQ0.uni = funcQ0.uni, funcQ0.biv = funcQ0.biv, funcQ0.tri = funcQ0.tri, funcQ0.quad = funcQ0.quad)
		class(out) <- "makeRandomData"

		# get summary
		bigObs <- remakeRandomData(n = 1e5, object = out)
		currR2 <- 1 - mean((bigObs$Y - bigObs$Q0)^2)/var(bigObs$Y)
		# cat("Current R2 = ", currR2)
	}

	return(out)
}
