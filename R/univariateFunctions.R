#' linUni
#'
#' Simulate a linear univariate relationship
#'
#' @param x A vector
#' @param coef A numeric (supplied by \code{linUniParm})
#'
#' @return \code{coef*x}
#' @export

linUni <- function(x, coef){
	# coef <- runif(1, -2, 2)
	return(coef*x)
}

#' linUniParm
#'
#' Draw random coefficient for linear univariate relationship
#' @export
linUniParm <- function(coefLower = -10, coefUpper = 10, ...){
	return(list(coef = runif(1, coefLower, coefUpper)))
}

#' linSplineUni
#' @export
linSplineUni <- function(x, nKnot, knotLoc, slopes){
	y <- rep(0, length(x))
	kl <- c(min(x),knotLoc)
	# make basis
	for(k in 1:(nKnot+1)){
		xb <- pmax(0, x - kl[k])
		y <- y + xb * slopes[k]
	}
	return(y)
}

#' linSplineUniParm
#' @export
linSplineUniParm <- function(x, slopeLower = -10, slopeUpper = 10,
                             nKnotLower = 1, nKnotUpper = 5,
                             quantLower = 0.05, quantUpper = 0.95){
	nKnot <- round(runif(1, nKnotLower - 0.5, nKnotUpper + 0.5))
	slopes <- runif(nKnot + 1, slopeLower, slopeUpper)
	knotLoc <- sort(quantile(x, p = runif(nKnot, quantLower, quantUpper)))
	return(list(nKnot = nKnot, slopes = slopes, knotLoc = knotLoc))
}


#' cubicSplineUni
#' @export
cubicSplineUni <- function(x, main, nKnot, knotLoc, slopes){
	kl <- c(min(x),knotLoc)
	y <- main[1]*x + main[2]*x^2 + main[3]*x^3
	# make basis
	for(k in 1:(nKnot+1)){
		xb <- pmax(0, x - kl[k])^3
		y <- y + xb * slopes[k]
	}
	return(y)
}

#' cubicSplineUniParm
#' @export
cubicSplineUniParm <- function(x, mainLower = -0.5, mainUpper = 0.5,
                              splineLower = -1, splineUpper = 1,
                             nKnotLower = 1, nKnotUpper = 20,
                             quantLower = 0.05, quantUpper = 0.95){
	main <- runif(3, mainLower, mainUpper)
	nKnot <- round(runif(1, nKnotLower - 0.5, nKnotUpper + 0.5))
	slopes <- runif(nKnot + 1, splineLower, splineUpper)
	knotLoc <- sort(quantile(x, p = runif(nKnot, quantLower, quantUpper)))
	return(list(main = main, nKnot = nKnot, slopes = slopes, knotLoc = knotLoc))
}


#' polyUni
#'
#' Simulate a polynomial univariate relationship
#'
#' @param x A vector
#' @param npoly A numeric indicating the number of polynomial terms to simulate
#' @param coef A vector indicating the coefficient for each polynomial term
#' @param deg A vector indicating the degree of each polynomial term
#' @param coef A numeric (supplied by \code{linUniParm})
#'
#' @export
#' @return \code{coef[1]*x^deg[1] + ... + coef[npoly]*x^deg[npoly]}

polyUni <- function(x, npoly, coef, deg){
	# npoly <- round(runif(1,0.5,3.5))
	# coef <- runif(npoly, -1, 1)
	# deg <- round(runif(npoly, 0.5, 3.5))
	tmp <- rep(0, length(x))
	for(i in 1:npoly){
		tmp <- tmp + coef[i]*x^deg[i]
	}
	return(tmp)
}

#' polyUniParm
#'
#' Draw random coefficient for polynomial relationship univariate relationship.
#' See source code for distributions.
#'
#' @export

polyUniParm <- function(npolyLower = 1, npolyUpper = 4, coefLower = -3, coefUpper = 3,
                        degLower = 1, degUpper = 4, ...){
	npoly <- sample(npolyLower:npolyUpper, 1)
	coef <- runif(npoly, coefLower, coefUpper)
	deg <- sample(1:4, npoly, replace = TRUE)
	return(list(npoly = npoly, coef = coef, deg = deg))
}


#' sinUni
#'
#' Simulate a sinusoidal univariate relationship
#'
#' @param x A vector
#' @param p The periodicity of the sin function
#' @param amp The amplitude of the sin function
#'
#' @return \code{amp*sin(p*x)}
#' @export
sinUni <- function(x, p, amp){
	w <- amp*sin(p*x)
	return(w)
}

#' sinUniParm
#'
#' Draw random periodicity and amplitude for sinusoidal relationship univariate relationship.
#' See source code for distributions.
#' @export
sinUniParm <- function(pLower = -4, pUpper = 4, ampLower = -4, ampUpper = 4, ...){
	p <- runif(1, pLower, pUpper)
	amp <- runif(1, ampLower, ampUpper)
	return(list(p = p, amp = amp))
}

#' jumpUni
#'
#' Simulate a jump function univariate relationship
#'
#' @param x A vector
#' @param njump The number of jumps the function takes
#' @param jumpLoc A vector of the locations of those jumps
#' @param jumpVal A vector of length njump indicating the value of the function at each jump
#'
#' @return The value of the jump function.
#' @export

jumpUni <- function(x, njump, jumpLoc, jumpVal){
	tmp <- rep(0, length(x))
	for(i in 1:(length(jumpLoc)-1)){
		ind <- (x <= jumpLoc[i+1] & x >= jumpLoc[i])
		tmp[ind] <- jumpVal[i]
	}
	return(tmp)
}

#' jumpUniParm
#'
#' Generate parameters for \code{jumpUni}. Here jumpLoc is a function that takes x and njump as input and outputs
#' the locations of jumps. By contrast, when jumpLoc is fed into jumpUni, it is just a vector
#' of numbers. See source code for how exactly it works.
#' @export

jumpUniParm <- function(x, njumpLower = 1, njumpUpper = 8,
                        jumpLower = -4, jumpUpper = 4,
                        jumpLoc = function(x,njump){
                        	quantile(x, seq(0,1,length=njump+1))
                        }){
	njump <- round(runif(1, njumpLower - 0.5, njumpUpper + 0.5))
	jumpL <- do.call("jumpLoc", args = list(x = x, njump = njump))
	jumps <- runif(njump, jumpLower, jumpUpper)
	return(list(
	   njump = njump, jumpLoc = jumpL, jumpVal = jumps
	))
}

#' qGammaUni
#'
#' Simulate a gamma cdf relationship with expit(x)
#' @param x A vector
#' @param coef The coefficient in front of the gamma cdf
#' @param a The shape parameter
#' @param b The rate parameter
#' @export
qGammaUni <- function(x, coef, a, b, ...){
	return(coef*qgamma(plogis(x), a, b))
}

#' qGammaUniParm
#'
#' Generate parameters for \code{qGammaUni}.
#' @export
#'
qGammaUniParm <- function(coefLower = -4, coefUpper = 4, aLower = 0.25, aUpper = 7,
                          bLower = 0.25, bUpper = 7,...){
	coef <- runif(1, coefLower, coefUpper)
	a <- runif(1, aLower, aUpper)
	b <- runif(1, bLower, bUpper)
	return(list(coef = coef, a = a, b = b))
}

#' dNormUni
#'
#' Simulate a normal pdf relationship with x
#' @param x A vector
#' @param coef The coefficient in front of the normal pdf
#' @param mu The mean of the normal pdf
#' @param sig The sd of the normal pdf
#'
#' @export
dNormUni <- function(x, coef, mu, sig){
	return(coef*dnorm(x, mu, sig))
}

#' dNormUniParm
#'
#' Generate parameters for \code{dNormUni}
#' @export
dNormUniParm <- function(...,coefLower = -4, coefUpper = 4,
                         muLower=-5, muUpper=5, sigLower=0.25, sigUpper=4){
	coef <- runif(1, coefLower, coefUpper)
	mu <- runif(1, muLower, muUpper)
	sig <- runif(1, sigLower, sigUpper)
	return(list(coef = coef, mu = mu, sig = sig))
}

#' pLogisUni
#'
#' Simulate an expit relationship with x
#' @param x A vector
#' @param coef The coefficient in front of the normal pdf
#' @param mult The multiplier inside the expit function
#' @export
pLogisUni <- function(x, coef, mult, loc, scale){
	return(coef*plogis(x*mult, loc, scale))
}

#' pLogisUniParm
#'
#' Parameters for \code{pLogisUni}
#' @export
#'
pLogisUniParm <- function(...,coefLower=-4, coefUpper=4,
                          multLower=-4, multUpper=4,
                          locLower = -4, locUpper = 4,
                          scaleLower = 0.25, scaleUpper = 4){
	coef <- runif(1, coefLower, coefUpper)
	mult <- runif(1, multLower, multUpper)
    loc <- runif(1, locLower, locUpper)
    scale <- runif(1, scaleLower, scaleUpper)
	return(list(coef = coef, mult = mult, loc = loc, scale = scale))
}

#' dNormMixUni
#'
#' Simulate a mixture of normal pdfs relationship with x
#' @param x A vector
#' @param coef1 Multiplier in front of first pdf
#' @param coef2 Multiplier in front of second pdf
#' @param mu1 Mean of first pdf
#' @param mu2 Mean of second pdf
#' @param sig1 SD of first pdf
#' @param sig2 SD of second pdf
#' @export
dNormMixUni <- function(x, coef1, coef2, mu1, mu2, sig1, sig2){
	return(coef1*dnorm(x, mu1, sig1) + coef2*dnorm(x, mu2, sig2))
}

#' dNormMixUniParm
#'
#' Generate parameters for \code{dNormMixUni}
#' @export
dNormMixUniParm <- function(...,coef1Lower = -4, coef1Upper = 4,
                            coef2Lower = -4, coef2Upper = 4,
                            mu1Lower = -10, mu1Upper = 10,
                            mu2Lower = -10, mu2Upper = 10,
                            sig1Lower = 0.5, sig1Upper = 4,
                            sig2Lower = 0.5, sig2Upper = 4
                            ){
	return(
       list(coef1 = runif(1, coef1Lower, coef1Upper),
            coef2 = runif(1, coef2Lower, coef2Upper),
            mu1 = runif(1, mu1Lower, mu1Upper),
            mu2 = runif(1, mu2Lower, mu2Upper),
			sig1 = runif(1,sig1Lower,sig1Upper),
           	sig2 = runif(1,sig2Lower,sig2Upper)
            )
   )
}
