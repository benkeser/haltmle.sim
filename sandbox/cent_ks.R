#! /usr/bin/env Rscript 
# This file was used to submit the simulation files to 
# a slurm-based Unix system. Using the sce.sh shell script
# one can submit each simulation in sequence. First, data files are
# created for each simulation. Those data files are then analyzed 
# in the 'run' execution. Then the results are collated in the 'merge'
# execution. 

# get environment variables
MYSCRATCH <- Sys.getenv('MYSCRATCH')
RESULTDIR <- Sys.getenv('RESULTDIR')
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
MYSCRATCH[is.na(MYSCRATCH)] <- '.'
RESULTDIR[is.na(RESULTDIR)] <- '.'
STEPSIZE[is.na(STEPSIZE)] <- 1
TASKID[is.na(TASKID)] <- 0

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

# # load packages
library(arm)
library(plyr)
library(gam, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2")
library(caret)
library(haltmle.sim, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2")
library(Rsolnp, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2")
library(future, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2")
library(cvma, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2")
library(hal9001, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2")
library(drtmle, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2")
library(SuperLearner)
# full parm
ns <- c(250, 500, 1000, 2000)
bigB <- 1000


# # # simulation parameters
parm <- expand.grid(seed=1:bigB,
                    n=ns)

# parm <- find_missing_files(tag = "ks",
#                            # parm needs to be in same order as 
#                            # file saves -- should make this more general...
#                            parm = c("n", "seed"),
#                            parm_values = list(seed = 1:1000,
#                                               n = c(200,1000,5000)))
# save(parm, file = "~/haltmle.sim/scratch/remain_ks_sims.RData")

# load("~/haltmle.sim/scratch/remain_ks_sims.RData")


# directories to save in 
saveDir <- "~/haltmle.sim/out/"
scratchDir <- "~/haltmle.sim/scratch/"

# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
  # for testing, useful to just run the first job
  # to make sure everything saves correctly
  # cat(1)
}

make_ks_mod <- function(n){
  z1 <- truncnorm::rtruncnorm(n, a = 0.5, b = 2.5, mean = 1.5)
  z2 <- truncnorm::rtruncnorm(n, a = -4, b = 4)
  z3 <- truncnorm::rtruncnorm(n, a = -4, b = 4)
  z4 <- truncnorm::rtruncnorm(n, a = -4, b = 4)

  w1 <- exp(z1/2)
  w2 <- z2/(1+exp(z1)) + 10
  w3 <- (z1*z3/25 + 0.6)^3
  w4 <- (z2 + z4 + 20)^2

  Y <- 210 + 27.4*z1 + 13.7*z2 + 13.7*z3 + 13.7*z4 + rnorm(n)
  g0 <- plogis(-z1 + 0.5*z2 - 0.25*z3 - 0.1*z4)
  A <- rbinom(n, 1, g0)
  return(list(Y = Y, A = A, W = data.frame(W1 = w1, W2 = w2, W3 = w3, W4 = w4),
              Z = data.frame(Z1 = z1, Z2 = z2, Z3 = z3, Z4 = z4)))
}

make_ks <- function(n){
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  z3 <- rnorm(n)
  z4 <- rnorm(n)

  w1 <- exp(z1/2)
  w2 <- z2/(1+exp(z1)) + 10
  w3 <- (z1*z3/25 + 0.6)^3
  w4 <- (z2 + z4 + 20)^2

  Y <- 210 + 27.4*z1 + 13.7*z2 + 13.7*z3 + 13.7*z4 + rnorm(n)
  g0 <- plogis(-z1 + 0.5*z2 - 0.25*z3 - 0.1*z4)
  A <- rbinom(n, 1, g0)
  return(list(Y = Y, A = A, W = data.frame(W1 = w1, W2 = w2, W3 = w3, W4 = w4),
              Z = data.frame(Z1 = z1, Z2 = z2, Z3 = z3, Z4 = z4)))
}

get_var_eif <- function(n = 1e6){
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  z3 <- rnorm(n)
  z4 <- rnorm(n)

  w1 <- exp(z1/2)
  w2 <- z2/(1+exp(z1)) + 10
  w3 <- (z1*z3/25 + 0.6)^3
  w4 <- (z2 + z4 + 20)^2
  Q0 <- 210 + 27.4*z1 + 13.7*z2 + 13.7*z3 + 13.7*z4
  Y <-  Q0 + rnorm(n)
  g0 <- plogis(-z1 + 0.5*z2 - 0.25*z3 - 0.1*z4)
  A <- rbinom(n, 1, g0)
  EIF <- (2*A - 1)/ifelse(A==1, g0, 1-g0) * (Y - Q0)
  return(var(EIF))
}

# execute prepare job ##################
if (args[1] == 'prepare') {

}

# execute parallel job #################################################
if (args[1] == 'run') {
  if (length(args) < 2) {
    stop("Not enough arguments. 'run' needs a second argument 'id'")
  }
  id <- as.numeric(args[2])
  print(paste(Sys.time(), "arrid:" , id, "TASKID:",
              TASKID, "STEPSIZE:", STEPSIZE))
  for (i in (id+TASKID):(id+TASKID+STEPSIZE-1)) {
    print(paste(Sys.time(), "i:" , i))

    # load parameters

    print(parm[i,])

    set.seed(parm$seed[i])
    dat <- make_ks(n=parm$n[i])

    algo <- c("SL.mean",
              "SL.hal9002",
              "SL.earth.cv",
              "SL.glm",
              "SL.bayesglm", 
              "SL.earth",
              "SL.step.interaction",
              "SL.gam", 
              "SL.dbarts.mod",
              "SL.gbm.caretMod",
              "SL.rf.caretMod",
              "SL.mean")
    # algo <- c("SL.glm","SL.bayesglm","SL.gam")
    # fit super learner with all algorithms
    set.seed(12314)
    out <- get_all_ates(Y = dat$Y, A = dat$A, W = dat$W, gtol = 0.025, 
                        V = 6, learners = algo, remove_learner = "SL.hal9002",
                        which_dr_tmle = c("full_sl", "cv_full_sl", "SL.hal9002", "cv_SL.hal9002")) #,
                        # which_dr_tmle = "full_sl")


    save(out, file=paste0(saveDir,"ksfinalmod_n=",parm$n[i],"_seed=",parm$seed[i],
                          ".RData"))
    }
}

# merge job ###########################
if (args[1] == 'merge') {
  # get_mean_rslt <- function(est_rslt){
  #   mean_rslt <- plyr::ddply(est_rslt, .(n), colMeans)
  #   return(mean_rslt)
  # }

  format_result <- function(out){
    tmle_os_list <- lapply(out, function(l){
      # browser()
      cv_est <- grepl("cv_", names(l))
      for(i in 1:length(l)){
        l[[i]] <- data.frame(as.list(l[[i]]))
      }
      for(i in 1:length(l)){
        if(cv_est[i]){
          l[[i]]$cv_ci_l <- l[[i]]$est - qnorm(0.975)*l[[i]]$se
          l[[i]]$cv_ci_u <- l[[i]]$est + qnorm(0.975)*l[[i]]$se
          l[[i]]$cov_cv_ci <- l[[i]]$cv_ci_l < truth & l[[i]]$cv_ci_u > truth
        }else{
          l[[i]]$ci_l <- l[[i]]$est - qnorm(0.975)*l[[i]]$se
          l[[i]]$ci_u <- l[[i]]$est + qnorm(0.975)*l[[i]]$se
          l[[i]]$cov_ci <- l[[i]]$ci_l < truth & l[[i]]$ci_u > truth
          l[[i]]$cv_ci_l <- l[[i]]$est - qnorm(0.975)*l[[i+1]]$se
          l[[i]]$cv_ci_u <- l[[i]]$est + qnorm(0.975)*l[[i+1]]$se
          l[[i]]$cov_cv_ci <- l[[i]]$cv_ci_l < truth & l[[i]]$cv_ci_u > truth
        }
      }
      return(l)
    })
    return(tmle_os_list)
  }
  truth <- 0
  all_files <- list.files("~/haltmle.sim/out")
  ks_files <- all_files[grepl("ksfinalmod",all_files)]
  logistic_tmle_rslt <- matrix(nrow = length(ks_files), ncol = 78 + 1)
  linear_tmle_rslt <- matrix(nrow = length(ks_files), ncol = 78 + 1)
  onestep_rslt <- matrix(nrow = length(ks_files), ncol = 78 + 1)
  drtmle_rslt <- matrix(nrow = length(ks_files), ncol = 26 + 1)
  for(i in seq_along(ks_files)){
    # get sample size
    # this_n <- as.numeric(strsplit(strsplit(ks_files[i], "_")[[1]][2], "n=")[[1]][2])
    this_n <- as.numeric(strsplit(strsplit(ks_files[i], "_")[[1]][3], "n=")[[1]][2])
    # load file
    load(paste0("~/haltmle.sim/out/",ks_files[i]))
      # format this file
    tmp <- format_result(out)
    # add results to rslt
    logistic_tmle_rslt[i,] <- c(this_n, unlist(tmp[[1]]))
    linear_tmle_rslt[i,] <- c(this_n,unlist(tmp[[2]]))
    onestep_rslt[i,] <- c(this_n,unlist(tmp[[3]]))
    drtmle_rslt[i,] <- c(this_n,unlist(tmp[[4]]))
  }
  col_names_1 <- c("n", names(unlist(tmp[[1]],use.names = TRUE)))
  col_names_2 <- c("n", names(unlist(tmp[[4]],use.names = TRUE)))
  logistic_tmle_rslt <- data.frame(logistic_tmle_rslt)
  colnames(logistic_tmle_rslt) <- col_names_1
  linear_tmle_rslt <- data.frame(linear_tmle_rslt)
  colnames(linear_tmle_rslt) <- col_names_1
  onestep_rslt <- data.frame(onestep_rslt)
  colnames(onestep_rslt) <- col_names_1
  drtmle_rslt <- data.frame(drtmle_rslt)
  colnames(drtmle_rslt) <- col_names_2
  rslt <- list(log_tmle = logistic_tmle_rslt,
               lin_tmle = linear_tmle_rslt,
               onestep = onestep_rslt,
               drtmle = drtmle_rslt)
  save(rslt, file = "~/haltmle.sim/out/allOut_finalmod.RData")
  
  # coverage
  by(rslt$lin_tmle, rslt$lin_tmle$n, function(x){
    colMeans(x[, grep("cov",colnames(x))
             # c("full_sl.cov_ci","full_sl.cov_cv_ci","cv_full_sl.cov_cv_ci",
             #      # "SL.rf.caretMod.cov_ci", "SL.glm.cov_ci",
             #      "SL.hal9001.cov_ci","SL.hal9001.cov_cv_ci","cv_SL.hal9001.cov_cv_ci")
             ])
  })  

  by(rslt$onestep, rslt$onestep$n, function(x){
    return(colMeans(x[, grep(".est",colnames(x))] * sqrt(x$n)))
  })



  est <- "SL.rf.caretMod"

  by(linear_tmle_rslt, logistic_tmle_rslt$n, function(x){
    colMeans(x[,paste0(c("","cv_"),est,".est")])
  })



  # by(logistic_tmle_rslt, logistic_tmle_rslt$n, function(x){
  #   c(sd(x$full_sl.est), mean(x$full_sl.se), sd(x$cv_full_sl.est), mean(x$cv_full_sl.se))
  # })




  # by(logistic_tmle_rslt, logistic_tmle_rslt$n, function(x){
  #   # colMeans(x[,c("full_sl.est","cv_full_sl.est",
  #   #               "SL.hal9001.est","cv_SL.hal9001.est"
  #   #               )])    
  #   colMeans(x[,grep(".est",colnames(x))])
  # })  


  # by(logistic_tmle_rslt, logistic_tmle_rslt$n, function(x){
  #   colMeans(x)
  # })

  # by(linear_tmle_rslt, linear_tmle_rslt$n, function(x){
  #   colMeans(x[,c("full_sl.cov_ci","full_sl.cov_cv_ci","cv_full_sl.cov_cv_ci",
  #                 "SL.hal9001.cov_ci","SL.hal9001.cov_cv_ci","cv_SL.hal9001.cov_cv_ci"
  #                 )])
  # })
  # by(onestep_rslt, onestep_rslt$n, function(x){
  #   colMeans(x[,c("full_sl.cov_ci","full_sl.cov_cv_ci","cv_full_sl.cov_cv_ci",
  #                 "SL.hal9001.cov_ci","SL.hal9001.cov_cv_ci","cv_SL.hal9001.cov_cv_ci"
  #                 )])
  # })
  # by(drtmle_rslt, drtmle_rslt$n, function(x){
  #   colMeans(x[,c("full_sl.cov_ci","full_sl.cov_cv_ci","cv_full_sl.cov_cv_ci",
  #                 "SL.hal9001.cov_ci","SL.hal9001.cov_cv_ci","cv_SL.hal9001.cov_cv_ci"
  #                 )])
  # })
}


