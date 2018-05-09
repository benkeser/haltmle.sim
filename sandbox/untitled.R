set.seed(1107)
makeData <- function(n, prob = 0.25, equal = TRUE){
  adjustVars <- data.frame(w1 = rnorm(n), w2 = rbinom(n,1,0.5))
  if (equal == T){
    trt <- sample(c(0:3), n, replace = T)
  }
  else {trt <- rbinom(n, 3, prob)}
  ftype <- rbinom(n, 1, 0.5) + 1 
  ctime <- 1 + rgeom(n, plogis(-4)) # random censoring
  ftime <- 1 + rgeom(n, plogis(0.5*adjustVars$w1 - adjustVars$w2*trt))
  time <- pmin(ctime, ftime)
  ftype[ctime < ftime] <- 0
  return(list(adjustVars = adjustVars, trt = trt, ftime = time, ftype = ftype))
}
dat <- makeData(1000)

debug(mean_tmle)
fit1 <- mean_tmle(ftime = dat$ftime, ftype = dat$ftype,
                  trt = dat$trt, adjustVars = dat$adjustVars,
                  ftypeOfInterest = 1,
                  glm.trt = "1", t0 = 6,
                  glm.ftime = "w1 + trt*w2",
                  glm.ctime = "w1 + trt + I(w1^2) + w2", msm.formula = "ftype + trt*w2",
                  msm.family = "binomial")


fit1