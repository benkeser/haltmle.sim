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

learners <- c("SL.hal9002",
              "SL.glm",
              "SL.bayesglm", 
              "SL.earth",
              "SL.step.interaction",
              "SL.gam", 
              "SL.dbarts.mod",
              "SL.gbm.caretMod",
              "SL.rf.caretMod",
              "SL.rpart.caretMod", 
              "SL.mean")

V <- 10
train_cols1 <- 1:V
train_cols2 <- combn(V, V-1)
train_cols3 <- combn(V, V-2)

tasks <- c("fit_or","fit_ps")

seeds <- 1:10000

nuisance_parm <- expand.grid(
                  learner = learners,
                  train_cols = c(split(train_cols2, col(train_cols2)),
                                 split(train_cols3), col(train_cols3),
                                 train_cols1),
                  seed = seeds, task = tasks)




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
    dat <- haltmle.sim:::makeRandomData(n=parm$n[i], maxD = 8,
                                        minObsA = 30,
                                        minR2 = parm$range_R2[[i]][1],
                                        maxR2 = parm$range_R2[[i]][2])    

    dat <- haltmle.sim:::makeRandomDataT(n=100, maxD = 2,
                                        # minObsA = 30,
                                        minR2 = 0.05,
                                        maxR2 = 0.95)

    save(dat, file=paste0(scratchDir,"draw_n=",parm$n[i],"_seed=",parm$seed[i],
                          "_r2=",parm$range_R2[[i]][1],".RData"))
    print("data saved")

    algo <- c("SL.hal9001",
              "SL.glm",
              "SL.bayesglm", 
              "SL.earth",
              "SL.step.interaction",
              "SL.gam", 
              "SL.dbarts.mod",
              "SL.gbm.caretMod",
              "SL.rf.caretMod",
              "SL.rpart.caretMod", 
              "SL.mean",
              "SL.kernelKnn")
    algo <- c("SL.glm","SL.hal9002")
    # fit super learner with all algorithms
    set.seed(parm$seed[i])
    dat$W <- data.frame(dat$W)
    colnames(dat$W) <- paste0("W",1:ncol(dat$W))
    debug(estimate_nuisance)
    debug(get_all_ates)
    debug(get_ctmle_fit)
    debug(predict_alllambda_SL.hal9002)
    out <- get_all_ates(Y = dat$Y$Y, A = dat$A$A, W = dat$W, 
                        V = 4, learners = algo, 
                        # remove_learner = "SL.hal9001",
                        which_Match = algo,
                        which_dr_tmle = algo,
                        which_dr_iptw = algo,
                        which_ctmle_g = "SL.hal9002")
    

    

    save(out, file=paste0(saveDir,"out_n=",parm$n[i],"_seed=",parm$seed[i],
                          "_r2=",parm$range_R2[[i]][1],".RData"))
    }
}

# merge job ###########################
if (args[1] == 'merge') {
  # format_result <- function(out){
  #   tmle_os_list <- lapply(out, function(l){
  #     # browser()
  #     cv_est <- grepl("cv_", names(l))
  #     for(i in 1:length(l)){
  #       l[[i]] <- data.frame(as.list(l[[i]]))
  #     }
  #     for(i in 1:length(l)){
  #       if(cv_est[i]){
  #         l[[i]]$cv_ci_l <- l[[i]]$est - qnorm(0.975)*l[[i]]$se
  #         l[[i]]$cv_ci_u <- l[[i]]$est + qnorm(0.975)*l[[i]]$se
  #         l[[i]]$cov_cv_ci <- l[[i]]$cv_ci_l < truth & l[[i]]$cv_ci_u > truth
  #       }else{
  #         l[[i]]$ci_l <- l[[i]]$est - qnorm(0.975)*l[[i]]$se
  #         l[[i]]$ci_u <- l[[i]]$est + qnorm(0.975)*l[[i]]$se
  #         l[[i]]$cov_ci <- l[[i]]$ci_l < truth & l[[i]]$ci_u > truth
  #         l[[i]]$cv_ci_l <- l[[i]]$est - qnorm(0.975)*l[[i+1]]$se
  #         l[[i]]$cv_ci_u <- l[[i]]$est + qnorm(0.975)*l[[i+1]]$se
  #         l[[i]]$cov_cv_ci <- l[[i]]$cv_ci_l < truth & l[[i]]$cv_ci_u > truth
  #       }
  #     }
  #     return(l)
  #   })
  #   return(tmle_os_list)
  # }
  # all_files <- list.files("~/haltmle.sim/out")
  # out_files <- all_files[grepl("out",all_files)][1:25]
  # data_summary <- matrix(nrow = length(out_files), ncol = 30 + 3)
  # logistic_tmle_rslt <- matrix(nrow = length(out_files), ncol = 182 + 3)
  # linear_tmle_rslt <- matrix(nrow = length(out_files), ncol = 182 + 3)
  # onestep_rslt <- matrix(nrow = length(out_files), ncol = 182 + 3)
  # drtmle_rslt <- matrix(nrow = length(out_files), ncol = 26 + 3)
  # for(i in seq_along(out_files)){
  #   # get sample size
  #   this_n <- as.numeric(strsplit(strsplit(out_files[i], "_")[[1]][2], "n=")[[1]][2])
  #   this_r2 <- as.numeric(strsplit(strsplit(strsplit(out_files[i], "_")[[1]][4], "r2=")[[1]][2],
  #                         ".RData")[[1]][1])
  #   this_seed <- as.numeric(strsplit(strsplit(out_files[i], "_")[[1]][3], "seed=")[[1]][2])
  #   # load file
  #   load(paste0("~/haltmle.sim/out/",out_files[i]))
  #   # load data set and compute summary
  #   load(paste0("~/haltmle.sim/scratch/",gsub("out","draw",out_files[i])))
  #   sum_dat <- summary(dat)[-c(13,14)]
  #   truth <- sum_dat$truth
  #   # format this file
  #   tmp <- format_result(out)
  #   # add results to rslt
  #   data_summary[i,] <- c(this_n, this_seed, this_r2, unlist(sum_dat)[2:31])
  #   logistic_tmle_rslt[i,] <- c(this_n, this_seed, this_r2, unlist(tmp[[1]]))
  #   linear_tmle_rslt[i,] <- c(this_n,this_seed, this_r2,unlist(tmp[[2]]))
  #   onestep_rslt[i,] <- c(this_n,this_seed, this_r2,unlist(tmp[[3]]))
  #   drtmle_rslt[i,] <- c(this_n,this_seed, this_r2,unlist(tmp[[4]]))
  # }
  # col_names_1 <- c("n","seed", "r2", names(unlist(tmp[[1]],use.names = TRUE)))
  # col_names_2 <- c("n","seed", "r2", names(unlist(tmp[[4]],use.names = TRUE)))
  # col_names_3 <- c("n","seed", "r2", names(unlist(sum_dat)[2:31]))
  # data_summary_rslt <- data.frame(apply(data_summary[,1:(ncol(data_summary)-1)], 2, as.numeric))
  # data_summary_rslt <- cbind(data_summary_rslt, data_summary[,ncol(data_summary)])
  # colnames(data_summary_rslt) <- col_names_3
  # logistic_tmle_rslt <- data.frame(logistic_tmle_rslt)
  # colnames(logistic_tmle_rslt) <- col_names_1
  # linear_tmle_rslt <- data.frame(linear_tmle_rslt)
  # colnames(linear_tmle_rslt) <- col_names_1
  # onestep_rslt <- data.frame(onestep_rslt)
  # colnames(onestep_rslt) <- col_names_1
  # drtmle_rslt <- data.frame(drtmle_rslt)
  # colnames(drtmle_rslt) <- col_names_2
  # rslt <- list(data_summary = data_summary_rslt,
  #              log_tmle = logistic_tmle_rslt,
  #              lin_tmle = linear_tmle_rslt,
  #              onestep = onestep_rslt,
  #              drtmle = drtmle_rslt)
  # save(rslt, file = "~/haltmle.sim/out/allOut_rdg.RData")

  load("~/Dropbox/R/haltmle.sim/results/allOut_ks.RData")
  # for computing absolute bias per unit information
  summary_all_bias <- function(rslt){
    tmp <- lapply(rslt[2:5], function(r){
      # find estimator columns
      est_col_idx <- grepl(".est", colnames(r))
      grbg <- sqrt(r$n) * abs(r[,est_col_idx] - rslt[[1]]$truth) / sqrt(rslt[[1]]$varEffIC)
      colMeans(grbg)
    })
    return(tmp)
  }

  # for computing coverage (only cv confidence intervals)
  summary_all_coverage <- function(rslt){
    tmp <- lapply(rslt[2:5], function(r){
      # find estimator columns
      ci_col_idx <- grepl(".cov_cv_ci", colnames(r))
      colMeans(r[,ci_col_idx])
    })
    return(tmp)
  }

  # for computing coverage (only cv confidence intervals)
  summary_all_mse <- function(rslt){
    tmp <- lapply(rslt[2:5], function(r){
      # find estimator columns
      col_idx <- grepl(".est", colnames(r))
      colMeans(r[,col_idx]^2)
    })
    return(tmp)
  }

  algo_frame <- function(){
    data.frame(est = c("full_sl","rm_sl",
              "SL.hal9002",
              "SL.earth.cv",
              "SL.glm",
              "SL.bayesglm", 
              # "SL.earth",
              "SL.step.interaction",
              # "SL.gam", 
              "SL.dbarts.mod",
              "SL.gbm.caretMod",
              "SL.rf.caretMod"),
              label = c("HAL-SL","SL","HAL","MARS","GLM","BGLM",
                        "STEPGLM","BART","GBM","RF"))
  }
  bias_correct_frame <- function(){
    data.frame(est = c("log_tmle", "lin_tmle","onestep","drtmle"),
               label = c("LOGTMLE", "LINTMLE", "OS", "DRTMLE"))
  }
  add_estimator_labels <- function(summary_measure_rslt){
    vec_smr <- unlist(summary_measure_rslt)
    cur_algo_names <- names(vec_smr)
    cur_bias_names <- names(summary_measure_rslt)
    algo_names <- rep(NA, length(unlist(summary_measure_rslt)))
    bias_correct_names <- rep(NA, length(unlist(summary_measure_rslt)))
    label_frame <- algo_frame()
    bias_frame <- bias_correct_frame()
    # get algorithm names
    for(i in seq_along(label_frame$est)){
      idx <- grep(label_frame$est[i], cur_algo_names)
      algo_names[idx] <- as.character(label_frame$label[i])
    }
    # get bias correct names
    for(i in seq_along(bias_frame$est)){
      idx <- grepl(bias_frame$est[i], cur_algo_names)
      idx_cv <- grepl("cv", cur_algo_names)
      bias_correct_names[idx & !idx_cv] <- as.character(bias_frame$label[i])
      bias_correct_names[idx & idx_cv] <- paste0("CV-",as.character(bias_frame$label[i]))
    }
    full_labels <- paste0(algo_names, " + ", bias_correct_names)
    return(full_labels)
  }
  make_ordered_plot <- function(rslt, # named list of results as above
                                stratify_conditions = NULL, # character of conditions to stratify on
                                summary_measure, # function to compute summary measure
                                comp_value = 0, # value to compare for ordering (0 for bias, 0.95 for coverage)
                                x_label = "", # label for x-axis
                                x_lim = c(0,1), ...# x limits
                                ){
    if(!is.null(stratify_conditions)){
      stratify_idx <- with(rslt$data_summary, which(eval(parse(text=stratify_conditions))))
      stratify_result <- lapply(rslt, function(r){
        r[stratify_idx,]
      })
    }else{
      stratify_result <- rslt
    }

    summary_measure_rslt <- do.call(summary_measure, args = list(rslt = stratify_result))
    estimator_labels <- add_estimator_labels(summary_measure_rslt)
    vec_smr <- unlist(summary_measure_rslt)
    order_vec_smr <- order(abs(vec_smr - comp_value))
    # plot
    par(oma = c(0,5.1,0,0), mar = c(5.1, 10.1, 2.1, 0.1))
    plot(x = 0, pch = "",
         y = 0,
         ylim = c(1,length(vec_smr)),
         yaxt = "n", bty = "n",
         xlab = x_label, ylab = "", 
         xlim = x_lim)

    library(RColorBrewer)
    log_col <- brewer.pal(n = 9, "Blues")[c(7,9)]
    lin_col <- brewer.pal(n = 9, "Greens")[c(7,9)]
    os_col <- brewer.pal(n = 9, "OrRd")[c(7,9)]
    drtmle_col <- brewer.pal(n = 9, "BuPu")[c(7,9)]
    ct <- length(estimator_labels) + 1
    for(i in seq_along(estimator_labels)){
      ct <- ct - 1
      cv_bool <- grepl("CV", estimator_labels[order_vec_smr][i])
      this_col <- if(grepl("LOGTMLE", estimator_labels[order_vec_smr][i])){
        log_col[1]
      }else if(grepl("LINTMLE", estimator_labels[order_vec_smr][i])){
        lin_col[1]  
      }else if(grepl("OS", estimator_labels[order_vec_smr][i])){
        os_col[1]
      }else if(grepl("DRTMLE", estimator_labels[order_vec_smr][i])){
        drtmle_col[1]
      }
      points(x = vec_smr[order_vec_smr][i], y = ct, pch = ifelse(cv_bool, 1, 4),
             col = this_col)
      axis(side = 2, lwd = 0, at = ct, label = estimator_labels[order_vec_smr][i],
           cex.axis = 0.25, las = 2, xpd = TRUE)
      abline(v = comp_value, lty = 3)
    }
  }

  new_rslt <- c(list(data_summary = data.frame(n = rslt[[1]]$n, truth = 0, varEffIC = 5.87)),
                rslt)

  make_ordered_plot(new_rslt, stratify_conditions = "n == 200", 
                    summary_measure = "summary_all_bias",
                    x_lim = c(0,100))  
  make_ordered_plot(new_rslt, stratify_conditions = "n == 1000", 
                    summary_measure = "summary_all_bias",
                    x_lim = c(0,100))  
  make_ordered_plot(new_rslt, stratify_conditions = "n == 5000", 
                    summary_measure = "summary_all_bias",
                    x_lim = c(0,100))

  make_ordered_plot(new_rslt, stratify_conditions = "n == 2000", 
                    summary_measure = "summary_all_mse",
                    x_lim = c(0,5), comp_value = 0)  
  plot_done()
  make_ordered_plot(new_rslt, stratify_conditions = "n == 1000", 
                    summary_measure = "summary_all_coverage",
                    x_lim = c(0,1), comp_value = 0.95)  
  make_ordered_plot(new_rslt, stratify_conditions = "n == 5000", 
                    summary_measure = "summary_all_coverage",
                    x_lim = c(0,1), comp_value = 0.95)



  #   # add color labels based on bias correction
  #   # lighter shade if CV, darker shade if non
  #   # 4 different colors

  # }
}


