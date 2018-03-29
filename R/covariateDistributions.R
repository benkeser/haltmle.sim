
#' uniformW
#'
#' Draw from uniform(lower, upper) distribution
#' @export
uniformW <- function(n, lower, upper,...){
	return(runif(n, lower, upper))
}
#' uniformWParm
#'
#' Generate parameters for \code{uniform}
#' @export

uniformWParm <- function(lowerLower = -10, lowerUpper = 0,
                        rangeLower = 0.5, rangeUpper = 10,...){
	lower <- runif(1, lowerLower, lowerUpper)
	upper <- lower + runif(1, rangeLower, rangeUpper)
	return(list(lower = lower, upper = upper))
}

#' bernoulliW
#'
#' Draw from Bernoulli(p)
#' @export
bernoulliW <- function(n,p,...){
	return(rbinom(n,1,p))
}
#' bernoulliParm
#'
#' Generate parameters for \code{bernoulli}
#' @export
bernoulliWParm <- function(pLower = 0.1, pUpper = 0.9,...){
	return(list(p = runif(1, pLower, pUpper)))
}

#' binomialW
#'
#' Draw from binomial(num, p)
#' @export
binomialW <- function(n,num,p,...){
	return(rbinom(n,num,p))
}
#' binomialWParm
#'
#' Generate parameters for \code{binomial}
#' @export
binomialWParm <- function(pLower = 0.1, pUpper = 0.9,
                          numLower = 2, numUpper = 20,...){
	return(list(p = runif(1, pLower, pUpper),
	            num = sample(numLower:numUpper, 1)))
}
#' normalW
#'
#' Draw from normal(mean,sd)
#' @export
normalW <- function(n,mean,sd,...){
	return(rnorm(n,mean,sd))
}
#' normalWParm
#'
#' Generate parameters for \code{normal}
#' @export

normalWParm <- function(meanLower = -2, meanUpper = 2,
                       sdLower = 0.25, sdUpper = 2,...){
	return(list(mean = runif(1, meanLower, meanUpper),
	            sd = runif(1, sdLower, sdUpper)))
}
#' gammaW
#'
#' Draw from gamma(a,b) distribution
#' @export
gammaW <- function(n,a,b,...){
	return(rgamma(n,a,b))
}
#' gammaWParm
#' Generate parameters for \code{gammaDist}
#' @export

gammaWParm <- function(aLower = 0.5, aUpper = 2.5,
                      bLower = 0.5, bUpper = 2.5,...){
	return(list(a = runif(1, aLower, aUpper),
	       b = runif(1, bLower, bUpper)))
}

#' normalWCor
#'
#' Draw normal variate with mean equal to sum of the selected variables
#' that it is correlated with
#' @export
normalWCor <- function(n,W=NULL,sd,whichCorr,...){
	if(!is.null(W)){
		return(rnorm(n,.rowSums(W),sd))
	}else{
		return(rnorm(n,0,sd))
	}
}
#' normalWCorParm
#'
#' Generate parameters for \code{normalCor} (smaller SD = more correlation)
#' @export
normalWCorParm <- function(W=NULL, sdLower = 0.25, sdUpper = 4,
                           nCorrLower = 1, nCorrUpper = ncol(W)){
	if(!is.null(W)){
		# randomly choose number of correlated variables
		nCorr <- round(runif(1, nCorrLower - 0.5, nCorrUpper + 0.5))
		# randomly choose which variables its correlated with
		whichCorr <- sample(1:ncol(W), nCorr)
	}else{
		whichCorr <- NA
	}
	return(list(sd = runif(1, sdLower, sdUpper),
	            whichCorr = whichCorr))
}

#' bernoulliWCor
#'
#' Draw bernoulli variable with prob x/r (to induce correlation)
#' @export
bernoulliWCor <- function(n,W = NULL,r,whichCorr,...){
	if(!is.null(W)){
		return(rbinom(n,1,plogis(.rowSums(W)/r)))
	}else{
		return(rbinom(n,1,plogis(0)))
	}
}

#' bernoulliWCorParm
#'
#' Generate parameters for \code{bernoulliCorParm} (smaller r = more correlation)
#' @export

bernoulliWCorParm <- function(W = NULL, rLower = -2, rUpper = 2,
                           nCorrLower = 1, nCorrUpper = ncol(W)){
	if(!is.null(W)){
		# randomly choose number of correlated variables
		nCorr <- sample(nCorrLower:nCorrUpper, 1)
		# randomly choose which variables its correlated with
		whichCorr <- sample(1:ncol(W), nCorr)
	}else{
		whichCorr <- NA
	}
	return(list(r = runif(1, rLower, rUpper),
	            whichCorr = whichCorr))
}

#' uniformWCor
#'
#' Draw variable that is x + uniform error (to induce correlation)
#' @export
uniformWCor <- function(n,W = NULL,lower,upper,whichCorr,...){
	if(!is.null(W)){
		return(.rowSums(W) + runif(n,lower,upper))
	}else{
		return(runif(n,lower,upper))
	}
}
#' uniformWParm
#' Generate parameters for \code{uniformParm} (smaller lower = more correlation)
#' @export

uniformWCorParm <- function(W = NULL, lowerLower = -5, lowerUpper = 0.25,
                        rangeLower = 0.25, rangeUpper = 10,
                        nCorrLower = 1, nCorrUpper = ncol(W)){
	if(!is.null(W)){
		# randomly choose number of correlated variables
		nCorr <- round(runif(1, nCorrLower - 0.5, nCorrUpper + 0.5))
		# randomly choose which variables its correlated with
		whichCorr <- sample(1:ncol(W), nCorr)
	}else{
		whichCorr <- NA
	}
	lower <- runif(1, lowerLower, lowerUpper)
	upper <- lower + runif(1, rangeLower, rangeUpper)
	return(list(lower = lower, upper = upper, whichCorr = whichCorr))
}

# distributions that are continuous but have point mass places
#' gammaPointMassW
#' @export
gammaPointMassW <- function(n,a,b,p,...){
	nPoint <- round(n*p)
	nGamma <- n - nPoint
	return(c(rgamma(nGamma, a, b), rep(0, nPoint)))
}

#' gammaWParm
#' Generate parameters for \code{gammaDist}
#' @export

gammaPointMassWParm <- function(aLower = 0.5, aUpper = 2.5,
                      bLower = 0.5, bUpper = 2.5,
                      pLower = 0.05, pUpper = 0.95, ...){
	return(list(a = runif(1, aLower, aUpper),
	       b = runif(1, bLower, bUpper),
	       p = runif(1, pLower, pUpper)))
}

#' normalPointMassW
#' @export
normalPointMassW <- function(n, mean, sd, p, ...){
	nPoint <- round(n*p)
	nNorm <- n - nPoint
	return(c(rnorm(nNorm, mean, sd), rep(0, nPoint)))
}

#' normalPointMassWParm
#'
#' Generate parameters for \code{normal}
#' @export
normalPointMassWParm <- function(meanLower = -2, meanUpper = 2,
                       sdLower = 0.25, sdUpper = 2,
                       pLower = 0.05, pUpper = 0.95, ...){
	return(list(mean = runif(1, meanLower, meanUpper),
	            sd = runif(1, sdLower, sdUpper),
	            p = runif(1, pLower, pUpper)))
}
# distributions that are binomial but have some fractional values
# in between integers
#' binomialFracW
#'
#' Draw from binomial(num, p)
#' @export
binomialFracW <- function(n,num,p,pFrac,...){
	nFrac <- round(n*pFrac)
	nBinom <- n - nFrac
	tmp1 <- rbinom(n,num,p)
	tmp2 <- c(runif(nFrac, 0, 1), rep(0, nBinom))
	return(tmp1 + tmp2)
}

#' binomialWParm
#'
#' Generate parameters for \code{binomial}
#' @export
binomialFracWParm <- function(pLower = 0.1, pUpper = 0.9,
                          numLower = 2, numUpper = 10,
                          pFracLower = 0.05, pFracUpper = 0.95, ...){
	return(list(p = runif(1, pLower, pUpper),
	            num = round(runif(1, numLower - 0.5, numUpper + 0.5)),
	            pFrac = runif(1, pFracLower, pFracUpper)))
}

#' .rowSums
#'
#' rowSums to handle one dimension properly
.rowSums <- function(x){
	if(dim(x)[2] == 1){
		unlist(x)
	}else{
		rowSums(x)
	}
}
