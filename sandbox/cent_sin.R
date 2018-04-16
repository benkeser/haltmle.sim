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
# ns <- c(200)
bigB <- 1000
four_period <- c(0.5, 5, 10)
two_period <- c(0.5, 10)

# # simulation parameters
parm_Q <- expand.grid(seed=1:bigB,
                    n=ns,
                    Q_period = four_period,
                    g_period = two_period)
parm_g <- expand.grid(seed=1:bigB,
                    n=ns,
                    Q_period = two_period,
                    g_period = four_period)
parm_g <- parm_g[-which(parm_g$g_period == 0.5 | parm_g$g_period == 10), ]

parm <- rbind(parm_Q, parm_g)
# colnames(parm)[3:4] <- c("Qp","gp")

# redo_parm <- find_missing_files(tag = "sin",
#                            # parm needs to be in same order as 
#                            # file saves -- should make this more general...
#                            parm = c("n", "seed", "Qp", "gp"),
#                            full_parm = parm)
# parm <- redo_parm
# save(parm, file = "~/haltmle.sim/scratch/remain_sin_sims.RData")
load("~/haltmle.sim/scratch/remain_sin_sims.RData")
names(parm)[3:4] <- c("Q_period","g_period")
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

make_sin <- function(n, Qp, gp){
  w1 <- runif(n, 0, 2*pi)
  # variation norm = 0.8*Qp
  Q0 <- 0.4*sin(Qp*w1) + 0.5 
  # variation norm = 0.8*gp
  g0 <- 0.4*sin(gp*w1) + 0.5

  A <- rbinom(n ,1, g0)
  Y <- rbinom(n, 1, Q0)

  return(list(W = data.frame(W1 = w1), A = A, Y = Y))
}

get_var_ic <- function(n = 1e6, Qp, gp){
  w1 <- runif(n, 0, 2*pi)
  # variation norm = 0.8*Qp
  Q0 <- 0.4*sin(Qp*w1) + 0.5 
  # variation norm = 0.8*gp
  g0 <- 0.4*sin(gp*w1) + 0.5

  A <- rbinom(n ,1, g0)
  Y <- rbinom(n, 1, Q0)
  IC <- (2*A - 1) / ifelse(A==1, g0, 1-g0) * (Y - Q0)
  return(var(IC))
}


# execute prepare job ##################
if (args[1] == 'prepare') {
  # for(i in 1:nrow(parm)){
  # }
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
    dat <- make_sin(n=parm$n[i], Qp = parm$Q_period[i], gp = parm$g_period[i])

    algo <- "SL.hal9002"
        
    # fit super learner with all algorithms
    out <- get_all_ates(Y = dat$Y, A = dat$A, W = dat$W, compute_superlearner = FALSE,
                        V = 6, learners = algo, remove_learner = NULL,
                        which_dr_tmle = c("SL.hal9002", "cv_SL.hal9002"))    

    save(out, file=paste0(saveDir,"sin_n=",parm$n[i],"_seed=",parm$seed[i],
                          "_Qp=",parm$Q_period[i],"_gp=",parm$g_period[i],".RData"))
    }
}

# merge job ###########################
if (args[1] == 'merge') {
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
  all_files <- list.files("~/haltmle.sim/out")
  sin_files <- all_files[grepl("sin",all_files)]
  logistic_tmle_rslt <- matrix(nrow = length(sin_files), ncol =13 + 4)
  linear_tmle_rslt <- matrix(nrow = length(sin_files), ncol =13 + 4)
  onestep_rslt <- matrix(nrow = length(sin_files), ncol =13 + 4)
  drtmle_rslt <- matrix(nrow = length(sin_files), ncol =13 + 4)
  for(i in seq_along(sin_files)){
    # get sample size
    this_n <- as.numeric(strsplit(strsplit(sin_files[i], "_")[[1]][2], "n=")[[1]][2])
    this_Qp <- as.numeric(strsplit(strsplit(sin_files[i], "_")[[1]][4], "Qp=")[[1]][2])
    this_gp <- as.numeric(strsplit(strsplit(strsplit(sin_files[i], "_")[[1]][5], "gp=")[[1]][2],
                          ".RData")[[1]][1])
    this_seed <- as.numeric(strsplit(strsplit(sin_files[i], "_")[[1]][3], "seed=")[[1]][2])
    # load file
    load(paste0("~/haltmle.sim/out/",sin_files[i]))
    # format this file
    tmp <- format_result(out)
    # add results to rslt
    logistic_tmle_rslt[i,] <- c(this_n, this_seed, this_Qp, this_gp, unlist(tmp[[1]]))
    linear_tmle_rslt[i,] <- c(this_n,this_seed, this_Qp, this_gp,unlist(tmp[[2]]))
    onestep_rslt[i,] <- c(this_n,this_seed, this_Qp, this_gp,unlist(tmp[[3]]))
    drtmle_rslt[i,] <- c(this_n,this_seed, this_Qp, this_gp,unlist(tmp[[4]]))
  }
  col_names_1 <- c("n","seed","Qp","gp", names(unlist(tmp[[1]],use.names = TRUE)))
  col_names_2 <- c("n","seed","Qp","gp", names(unlist(tmp[[4]],use.names = TRUE)))
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
  save(rslt, file = "~/haltmle.sim/out/allOut_sin.RData")
}

if(FALSE){
# locally
  load("~/Dropbox/R/haltmle.sim/results/allOut_sin.RData")  

  get_mean_rslt <- function(est_rslt){
    mean_rslt <- ddply(est_rslt, .(n, Qp, gp), colMeans)
    return(mean_rslt)
  }

  
  make_plot <- function(mean_rslt, col_name, xvar = "Qp"){
    cov_plot <- grepl("cov", col_name)
    if(cov_plot){
      yl <- c(0.5, 1)
    }else{
      yl <- range(mean_rslt[,col_name])
    }

    plot(0, 0, xlim = c(0, 10), ylim = yl, pch = "")
    stratvar <- ifelse(xvar == "Qp", "gp", "Qp")

    ct <- 0
    for(n in unique(mean_rslt$n)){
      ct <- ct + 1
      yvar_low <- mean_rslt[mean_rslt$n == n & mean_rslt[,stratvar] == 10, col_name]
      xvar_low <- mean_rslt[mean_rslt$n == n & mean_rslt[,stratvar] == 10, xvar]
      points(y = yvar_low, x = xvar_low, pch = ct)
    }
    abline(h = ifelse(cov_plot, 0.95, 0), lty = 3)
  }

  library(RColorBrewer)
  make_density_plot <- function(est_rslt, col_name, xvar = "Qp", est_lab = "",
                                dens_col = brewer.pal(n = 9, name = "Greens")[c(3,4,5,7)]){
    layout(matrix(c(1:6),nrow = 2, byrow = TRUE))
    par(mar = c(1.1,3.1,3.1,1.1), mgp = c(2, 0.5, 0), oma = c(0,2.5,2.5,0))
    stratvar <- ifelse(xvar == "Qp", "gp", "Qp")
    for(n in sort(unique(est_rslt$n))){
      vic <- get_var_ic(Qp = 0.5, gp = 0.5)
      yvar_low <- est_rslt[est_rslt$n == n & est_rslt[,stratvar] == 0.5 & est_rslt[,xvar] == 0.5, col_name]
      ct <- 1
      plot(density(sqrt(n)*yvar_low/sqrt(vic)), col = dens_col[ct], lwd = 2, ylim = c(0,0.45),
           main = paste0("n = ",n), xlab = "", xlim = c(-4,4), bty = "n")
      if(n == min(unique(est_rslt$n))){
        if(xvar == "Qp"){
          tit <- expression("||"*bar(Q)[0]*"||"[v])
        }else{
          tit <- expression("||"*g[0]*"||"[v])
        }
        legend(x = "topleft", bty = "n", col = c(dens_col, 1), lwd = 2, lty = c(1,1,1,1,3),
               legend = c(0.8, 1.6, 8, 16, "N(0,1)"), title = tit)
      }
      for(i in c(1, 5, 10)){
        if(xvar == "Qp"){
          vic <- get_var_ic(Qp = i, gp = 0.5)
        }else{
          vic <- get_var_ic(Qp = 0.5, gp = i)
        }
        ct <- ct + 1
        yvar_low <- est_rslt[est_rslt$n == n & est_rslt[,stratvar] == 0.5 & est_rslt[,xvar] == i, col_name]
        lines(density(sqrt(n)*yvar_low/sqrt(vic)), col = dens_col[ct], lwd = 2)
      }
      lines(y = dnorm(seq(-4,4,length=1e3)), x = seq(-4,4,length=1e3), col = 1, lwd = 2, lty = 3)
    }
    par(mar = c(3.1,3.1,1.1,1.1), mgp = c(2, 0.5, 0))
    for(n in sort(unique(est_rslt$n))){
      if(xvar == "Qp"){
        vic <- get_var_ic(Qp = 0.5, gp = 10)
      }else{
        vic <- get_var_ic(Qp = 10, gp = 0.5)        
      }
      yvar_low <- est_rslt[est_rslt$n == n & est_rslt[,stratvar] == 10 & est_rslt[,xvar] == 0.5, col_name]
      ct <- 1
      plot(density(sqrt(n)*yvar_low/sqrt(vic)), col = dens_col[ct], lwd = 2, ylim = c(0,0.45),
           main = "", 
           # xlab = expression(frac(psi[n] - psi[0], sqrt("Var(EIF) / n"))), 
           xlab = expression(sqrt(n)*"("*psi[n] - psi[0]*") / "*sqrt("var(EIF)")),
           xlim = c(-4,4),
           bty = "n")

      for(i in c(1, 5, 10)){
        if(xvar == "Qp"){
          vic <- get_var_ic(Qp = i, gp = 10)
        }else{
          vic <- get_var_ic(Qp = 10, gp = i)          
        }
        ct <- ct + 1
        yvar_low <- est_rslt[est_rslt$n == n & est_rslt[,stratvar] == 10 & est_rslt[,xvar] == i, col_name]
        lines(density(sqrt(n)*yvar_low/sqrt(vic)), col = dens_col[ct], lwd = 2)
      }
      lines(y = dnorm(seq(-4,4,length=1e3)), x = seq(-4,4,length=1e3), col = 1, lwd = 2, lty = 3)
    }
    if(xvar == "Qp"){
      mtext(outer = TRUE, side = 2, expression("||"*g[0]*"||"[v]*" 0.8"), at = 0.725, cex = 0.75, line = 0.9)
      mtext(outer = TRUE, side = 2, expression("||"*g[0]*"||"[v]*" = 16"), at = 0.275, cex = 0.75, line = 0.9)
    }else{
      mtext(outer = TRUE, side = 2, expression("||"*bar(Q)[0]*"||"[v]*" 0.8"), at = 0.725, cex = 0.75, line = 0.9)
      mtext(outer = TRUE, side = 2, expression("||"*bar(Q)[0]*"||"[v]*" = 16"), at = 0.275, cex = 0.75, line = 0.9)
    }
    mtext(outer = TRUE, side = 3, est_lab)
  }

  # sinusoidal simulation results
  pdf("~/Dropbox/R/haltmle.sim/results/sin_density_results.pdf", height = 4.2, width = 8.2)
  make_density_plot(rslt[[1]], col_name = "SL.hal9001.est", xvar = "Qp", est_lab = "HAL + logistic TMLE")
  make_density_plot(rslt[[1]], col_name = "SL.hal9001.est", xvar = "gp", est_lab = "HAL + logistic TMLE")
  make_density_plot(rslt[[2]], col_name = "SL.hal9001.est", xvar = "Qp", est_lab = "HAL + linear TMLE")
  make_density_plot(rslt[[2]], col_name = "SL.hal9001.est", xvar = "gp", est_lab = "HAL + linear TMLE")
  make_density_plot(rslt[[3]], col_name = "SL.hal9001.est", xvar = "Qp", est_lab = "HAL + one step")
  make_density_plot(rslt[[3]], col_name = "SL.hal9001.est", xvar = "gp", est_lab = "HAL + one step")
  make_density_plot(rslt[[4]], col_name = "SL.hal9001.est", xvar = "Qp", est_lab = "HAL + DR TMLE")
  make_density_plot(rslt[[4]], col_name = "SL.hal9001.est", xvar = "gp", est_lab = "HAL + DR TMLE")  
  make_density_plot(rslt[[1]], col_name = "cv_SL.hal9001.est", xvar = "Qp", est_lab = "HAL + logistic CV-TMLE")
  make_density_plot(rslt[[1]], col_name = "cv_SL.hal9001.est", xvar = "gp", est_lab = "HAL + logistic CV-TMLE")
  make_density_plot(rslt[[2]], col_name = "cv_SL.hal9001.est", xvar = "Qp", est_lab = "HAL + linear CV-TMLE")
  make_density_plot(rslt[[2]], col_name = "cv_SL.hal9001.est", xvar = "gp", est_lab = "HAL + linear CV-TMLE")
  make_density_plot(rslt[[3]], col_name = "cv_SL.hal9001.est", xvar = "Qp", est_lab = "HAL + CV-one step")
  make_density_plot(rslt[[3]], col_name = "cv_SL.hal9001.est", xvar = "gp", est_lab = "HAL + CV-one step")
  make_density_plot(rslt[[4]], col_name = "cv_SL.hal9001.est", xvar = "Qp", est_lab = "HAL + CV-DR TMLE")
  make_density_plot(rslt[[4]], col_name = "cv_SL.hal9001.est", xvar = "gp", est_lab = "HAL + CV-DR TMLE")
  dev.off()


  make_coverage_plot <- function(est_rslt, col_name, xvar = "Qp", est_lab = "",
                              dens_col = brewer.pal(n = 9, name = "Greens")[c(3,4,5,7)]){
  layout(matrix(c(1:6),nrow = 2, byrow = TRUE))
  par(mar = c(1.1,3.1,3.1,1.1), mgp = c(2, 0.5, 0), oma = c(0,2.5,2.5,0))
  stratvar <- ifelse(xvar == "Qp", "gp", "Qp")
  for(n in sort(unique(est_rslt$n))){
    # vic <- get_var_ic(Qp = 0.5, gp = 0.5)
    yvar_low <- est_rslt[est_rslt$n == n & est_rslt[,stratvar] == 0.5 & est_rslt[,xvar] == 0.5, col_name]
    ct <- 1
    plot(x = 1, y = mean(yvar_low), col = dens_col[ct], lwd = 2, ylim = c(0.5, 1),
         main = paste0("n = ",n), xlab = "", xlim = c(1,4), bty = "n", xaxt = "n",
         pch = 19, ylab = "Coverage probability (95% MC CI)")
    val <- mean(yvar_low)
    ci <- rep(val,2) + c(-1.96, 1.96)*sqrt(val*(1-val)/1000)
    segments(x0 = 1, y0 = ci[1], y1 = ci[2], col = dens_col[ct], lwd = 2)
    axis(side = 1, at = 1:4, labels = c(0.8, 1.6, 8, 16))
    if(n == min(unique(est_rslt$n))){
      if(xvar == "Qp"){
        tit <- expression("||"*bar(Q)[0]*"||"[v])
      }else{
        tit <- expression("||"*g[0]*"||"[v])
      }
      # legend(x = "topleft", bty = "n", col = c(dens_col, 1), lwd = 2, lty = c(1,1,1,1,3),
      #        legend = c(0.8, 1.6, 8, 16, "N(0,1)"), title = tit)
    }
    for(i in c(1, 5, 10)){
      ct <- ct + 1
      yvar_low <- est_rslt[est_rslt$n == n & est_rslt[,stratvar] == 0.5 & est_rslt[,xvar] == i, col_name]
      val <- mean(yvar_low)
      points(x = ct, y = val, col = dens_col[ct], pch = 19)
      ci <- rep(val,2) + c(-1.96, 1.96)*sqrt(val*(1-val)/1000)
      segments(x0 = ct, y0 = ci[1], y1 = ci[2], col = dens_col[ct], lwd = 2)
    }
    abline(h = 0.95, lty = 3)
  }
  par(mar = c(3.1,3.1,1.1,1.1), mgp = c(2, 0.5, 0))
  for(n in sort(unique(est_rslt$n))){
    yvar_low <- est_rslt[est_rslt$n == n & est_rslt[,stratvar] == 10 & est_rslt[,xvar] == 0.5, col_name]
    ct <- 1
    if(n == min(unique(est_rslt$n))){
      if(xvar == "Qp"){
        tit <- expression("||"*bar(Q)[0]*"||"[v])
      }else{
        tit <- expression("||"*g[0]*"||"[v])
      }
      # legend(x = "topleft", bty = "n", col = c(dens_col, 1), lwd = 2, lty = c(1,1,1,1,3),
      #        legend = c(0.8, 1.6, 8, 16, "N(0,1)"), title = tit)
    }
    val <- mean(yvar_low) 
    plot(x = 1, y = val, col = dens_col[ct], lwd = 2, ylim = c(0.5, 1),
         main = "", xlab = tit, xlim = c(1,4), bty = "n", xaxt = "n",
         pch = 19, ylab = "Coverage probability (95% MC CI)")
    ci <- rep(val,2) + c(-1.96, 1.96)*sqrt(val*(1-val)/1000)
    segments(x0 = 1, y0 = ci[1], y1 = ci[2], col = dens_col[ct], lwd = 2)

    axis(side = 1, at = 1:4, labels = c(0.8, 1.6, 8, 16))
    for(i in c(1, 5, 10)){
      ct <- ct + 1
      yvar_low <- est_rslt[est_rslt$n == n & est_rslt[,stratvar] == 10 & est_rslt[,xvar] == i, col_name]
      val <- mean(yvar_low)
      points(x = ct, y = mean(val), col = dens_col[ct], pch = 19)
      ci <- rep(val,2) + c(-1.96, 1.96)*sqrt(val*(1-val)/1000)
      segments(x0 = ct, y0 = ci[1], y1 = ci[2], col = dens_col[ct], lwd = 2)
    }
    abline(h = 0.95, lty = 3)
  }
  if(xvar == "Qp"){
    mtext(outer = TRUE, side = 2, expression("||"*g[0]*"||"[v]*" 0.8"), at = 0.725, cex = 0.75, line = 0.9)
    mtext(outer = TRUE, side = 2, expression("||"*g[0]*"||"[v]*" = 16"), at = 0.275, cex = 0.75, line = 0.9)
  }else{
    mtext(outer = TRUE, side = 2, expression("||"*bar(Q)[0]*"||"[v]*" 0.8"), at = 0.725, cex = 0.75, line = 0.9)
    mtext(outer = TRUE, side = 2, expression("||"*bar(Q)[0]*"||"[v]*" = 16"), at = 0.275, cex = 0.75, line = 0.9)
  }
  mtext(outer = TRUE, side = 3, est_lab)
  }

  pdf("~/Dropbox/R/haltmle.sim/results/sin_coverage_results.pdf", height = 4.2, width = 8.2)
  make_coverage_plot(rslt[[1]], col_name = "SL.hal9001.cov_cv_ci", xvar = "Qp", est_lab = "HAL + logistic TMLE")
  make_coverage_plot(rslt[[1]], col_name = "SL.hal9001.cov_cv_ci", xvar = "gp", est_lab = "HAL + logistic TMLE")
  make_coverage_plot(rslt[[2]], col_name = "SL.hal9001.cov_cv_ci", xvar = "Qp", est_lab = "HAL + linear TMLE")
  make_coverage_plot(rslt[[2]], col_name = "SL.hal9001.cov_cv_ci", xvar = "gp", est_lab = "HAL + linear TMLE")
  make_coverage_plot(rslt[[3]], col_name = "SL.hal9001.cov_cv_ci", xvar = "Qp", est_lab = "HAL + one step")
  make_coverage_plot(rslt[[3]], col_name = "SL.hal9001.cov_cv_ci", xvar = "gp", est_lab = "HAL + one step")
  make_coverage_plot(rslt[[4]], col_name = "SL.hal9001.cov_cv_ci", xvar = "Qp", est_lab = "HAL + DR TMLE")
  make_coverage_plot(rslt[[4]], col_name = "SL.hal9001.cov_cv_ci", xvar = "gp", est_lab = "HAL + DR TMLE")  
  make_coverage_plot(rslt[[1]], col_name = "cv_SL.hal9001.cov_cv_ci", xvar = "Qp", est_lab = "HAL + logistic CV-TMLE")
  make_coverage_plot(rslt[[1]], col_name = "cv_SL.hal9001.cov_cv_ci", xvar = "gp", est_lab = "HAL + logistic CV-TMLE")
  make_coverage_plot(rslt[[2]], col_name = "cv_SL.hal9001.cov_cv_ci", xvar = "Qp", est_lab = "HAL + linear CV-TMLE")
  make_coverage_plot(rslt[[2]], col_name = "cv_SL.hal9001.cov_cv_ci", xvar = "gp", est_lab = "HAL + linear CV-TMLE")
  make_coverage_plot(rslt[[3]], col_name = "cv_SL.hal9001.cov_cv_ci", xvar = "Qp", est_lab = "HAL + CV-one step")
  make_coverage_plot(rslt[[3]], col_name = "cv_SL.hal9001.cov_cv_ci", xvar = "gp", est_lab = "HAL + CV-one step")
  make_coverage_plot(rslt[[4]], col_name = "cv_SL.hal9001.cov_cv_ci", xvar = "Qp", est_lab = "HAL + CV-DR TMLE")
  make_coverage_plot(rslt[[4]], col_name = "cv_SL.hal9001.cov_cv_ci", xvar = "gp", est_lab = "HAL + CV-DR TMLE")
  dev.off()

}


