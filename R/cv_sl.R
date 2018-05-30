get_rm_sl_input <- function(split_Y, Yname, valid_folds, 
                         V, all_fit_tasks, all_fits, learners,
                         remove_index){
    # gives all splits used to compute this super learner
    if(!is.null(valid_folds)){
        sl_training_folds <- (1:V)[-valid_folds] 
    }else{
        # evaluates when getting outer super learner
        sl_training_folds <- 1:V
    }
    # gives a matrix whose columns correspond to the folds used 
    # to train each of the candidate 
    learner_training_folds <- combn(sl_training_folds, V - (length(valid_folds) + 1))

    # find indexes corresponding to these training sets
    # each column corresponds to the indexes where the corresponding
    # column of train_matrix were in the training folds
    # each row corresponds to a different candidate learner 
    # TO DO: better book keeping could allow this to be done without
    #        using search_fits
    fit_idx_list <- lapply(split(learner_training_folds, col(learner_training_folds)), 
                             cvma:::search_fits_for_training_folds,
                             fits = all_fit_tasks, y = Yname)
    fit_idx_list <- lapply(fit_idx_list, function(i){ i[-remove_index] })

    # let's assume that the method for obtaining weights will want an input
    # of Y, Y_hat (matrix from learners), and fold

    # need to put together matrix of predictions from learners
    Y_out <- lapply(split(learner_training_folds, col(learner_training_folds)), 
                    cvma:::get_Y_out, split_Y = split_Y, 
                    training_folds = sl_training_folds)

    # here x does get input to get_fold_out
    # valid_folds = folds in training_folds not used by the current fit
    fold_out <- lapply(split(learner_training_folds, col(learner_training_folds)), 
                       cvma:::get_fold_out, split_Y = split_Y, 
                       training_folds = sl_training_folds)

    pred_out <- mapply(learner_folds = split(learner_training_folds, col(learner_training_folds)), 
                       fit_idx = fit_idx_list,
                       cvma:::get_pred_out, 
                       MoreArgs = list(sl_folds = sl_training_folds,
                       all_fits = all_fits, 
                       learners = learners),
                       SIMPLIFY = FALSE)
    # reorganize into proper format
    out <- do.call(Map, c(list, list(fold_out, Y_out, pred_out)))
    out <- lapply(out, function(o){ names(o) <- c("valid_folds", "Y", "pred"); o})

    return(out)
}


get_rm_sl <- function(task, Y, V, all_fit_tasks, all_fits, folds, sl_control, learners,
                      remove_index){
    split_Y <- split(Y, folds)
    if(!all(1:V %in% task$training_folds)){
        valid_folds <- (1:V)[-task$training_folds]        
    }else{
        valid_folds <- NULL
    }
    # get input needed to compute sl ensemble weights 
    input <- get_rm_sl_input(split_Y = split_Y, valid_folds = valid_folds,
                          Yname = task$Yname, V = V, all_fit_tasks = all_fit_tasks, 
                          all_fits = all_fits, learners = learners, 
                          remove_index = remove_index)

    # get sl ensemble weights
    sl_weight <- do.call(sl_control$weight_fn, 
                          args = list(input = input, sl_control = sl_control))

    out <- list(training_folds = task$training_folds, Yname = task$Yname, sl_weight = sl_weight$weight)
    # get sl predictions on valid_folds folds by searching for when
    # (1:V)[-valid_folds] is training sample
    # sl_pred <- get_sl_pred(valid_folds = valid_folds, V = V, all_fit_tasks = all_fit_tasks,
    #                        all_fits = all_fits, y = Y_name, sl_weight = sl_weight,
    #                        split_Y = split_Y)
    return(out)
}

get_ate_cv_Q_pred <- function(Y, V, all_fit_tasks, all_fits, all_sl, folds, 
                             sl_control, return_learner_fits = TRUE,
                             learners, remove_index = NULL, compute_superlearner = TRUE){
	# 
  n <- length(Y)
	train_matrix <- combn(V, V-1)
	all_out <- lapply(split(train_matrix,col(train_matrix)), function(tr){
		learner_idx <- cvma:::search_fits_for_training_folds(fits = all_fit_tasks, 
                                          y = "Y", 
                                          training_folds = tr)
		if(compute_superlearner){
      if(!is.null(remove_index)){
  			learner_idx <- learner_idx[-remove_index]
  		}
  		sl_idx <- cvma:::search_fits_for_training_folds(fits = all_sl, 
                                            y = "Y", 
                                            training_folds = tr)
    }
		all_pred_setA <- Reduce(cbind,sapply(learner_idx, function(i){
			unlist(all_fits[[i]]$pred_setA)
		}, simplify = FALSE))

    if(compute_superlearner){
  		sl_pred_setA <- do.call(sl_control$ensemble_fn, args = list(pred = all_pred_setA,
  		                                                            weight = all_sl[[sl_idx]]$sl_weight))
    }else{
      sl_pred_setA <- NULL
    }
		return(list(learner_pred = all_pred_setA,
		            sl_pred = sl_pred_setA))
	})
  	idx <- unlist(split(1:n, folds)[V:1], use.names = FALSE)
  	cv_learner_pred <- matrix(NA, nrow = 2*length(Y), 
  	                          ncol = ifelse(is.null(remove_index), length(learners),
  	                                        length(learners) - length(remove_index)))
    tmp <- lapply(all_out, "[[", "learner_pred") # list 2x 
    cv_learner_pred_1 <- cv_learner_pred_0 <- matrix(NA, nrow = length(Y), ncol = ifelse(is.null(remove_index), length(learners),
                                            length(learners) - length(remove_index)))   
	if(length(learners) > 1){
    cv_pred_0 <- lapply(tmp, function(x){ nr <- nrow(x)/2; return(x[1:nr,]) })
    cv_pred_1 <- lapply(tmp, function(x){ nr <- nrow(x)/2; return(x[(nr+1):(2*nr),]) })
    cv_learner_pred_0[idx,] <- Reduce(rbind, cv_pred_0)
    cv_learner_pred_1[idx,] <- Reduce(rbind, cv_pred_1)
    cv_learner_pred <- rbind(cv_learner_pred_0, cv_learner_pred_1)
    # cv_learner_pred[c(rbind(idx, n + idx)),] <- Reduce(rbind, c(cv_pred_0, cv_pred_1))
  }else{
    cv_pred_0 <- lapply(tmp, function(x){ nr <- length(x)/2; return(x[1:nr]) })
    cv_pred_1 <- lapply(tmp, function(x){ nr <- length(x)/2; return(x[(nr+1):(2*nr)]) })
    cv_learner_pred_0[idx] <- Reduce(c, cv_pred_0)
    cv_learner_pred_1[idx] <- Reduce(c, cv_pred_1)
    cv_learner_pred <- matrix(c(cv_learner_pred_0, cv_learner_pred_1))
  }
	if(compute_superlearner){
    cv_sl_pred_0 <- cv_sl_pred_1 <- rep(NA, length(Y))
    tmp <- lapply(all_out, "[[", "sl_pred") # list 2x 
    cv_sl_pred_0[idx] <- unlist(lapply(tmp, function(x){ nr <- length(x)/2; return(x[1:nr]) }), use.names = FALSE)
    cv_sl_pred_1[idx] <- unlist(lapply(tmp, function(x){ nr <- length(x)/2; return(x[(nr+1):(2*nr)]) }), use.names = FALSE)
    cv_sl_pred <- c(cv_sl_pred_0, cv_sl_pred_1)
  }else{
    cv_sl_pred <- NULL
  }
    return(list(cv_learner_pred = cv_learner_pred,
                cv_sl_pred = cv_sl_pred))
}


#' Fit a learner on training folds and get predictions on validation folds
#' 
#' @param task A named list identifying what training_folds to fit the learner
#' on. The function returns predictions from this fit on the remaining folds (i.e.,
#' the validation folds). 
#' @param folds Vector identifying which fold observations fall into. 
#' @param Y A matrix or data.frame of outcomes.
#' @param X A matrix or data.frame of predictors.
#' @param sl_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' 
#' @return A named list with task and output of super learner wrapper fit. 

get_or_fit <- function(task, folds, W, A, Y, sl_control){
    train_idx <- folds %in% task$training_folds
    valid_idx <- !train_idx
    X <- data.frame(A = A, W)
    train_Y <- Y[train_idx]
    train_X <- X[train_idx, , drop = FALSE]
    if(sum(valid_idx) > 0){
      valid_setX <- rbind(data.frame(A = 0, W[valid_idx, , drop = FALSE]),
                       data.frame(A = 1, W[valid_idx, , drop = FALSE]))
      valid_X <- X[valid_idx, , drop = FALSE]
      n_valid_X <- sum(valid_idx)
    }else{
      # if learner being fit on all data, then return
      # predictions on training sample
      valid_setX <- rbind(data.frame(A = 0, W),
                       	  data.frame(A = 1, W))
      valid_X <- train_X  
      n_valid_X <- length(W[,1])
  	}
    # fit super learner wrapper
    fit <- do.call(task$SL_wrap, list(Y = train_Y, X = train_X,
                                 newX = rbind(valid_X, valid_setX),
                                 obsWeights = rep(1, length(train_Y)),
                                 family = sl_control$family))
    # retrieve predictions from observed A
    pred_obsA <- fit$pred[1:n_valid_X]
    pred_setA <- fit$pred[(n_valid_X+1):length(fit$pred)]
    # split up validation predictions
    if(sum(valid_idx) > 0){
      fit$pred <- split(pred_obsA, folds[valid_idx])
      fit$pred_setA <- split(pred_setA, folds[valid_idx])
    }else{
      fit$pred <- pred_obsA
      fit$pred_setA <- pred_setA
    }
    # return(fit)
    return(c(task, fit))
}


get_ps_fit <- function(task, folds, W, A, sl_control){
    train_idx <- folds %in% task$training_folds
    valid_idx <- !train_idx
    X <- W
    train_Y <- A[train_idx]
    train_X <- X[train_idx, , drop = FALSE]
    if(sum(valid_idx) > 0){
      valid_X <- X[valid_idx, , drop = FALSE]
      n_valid_X <- sum(valid_idx)
    }else{
      # if learner being fit on all data, then return
      # predictions on training sample
      valid_X <- train_X  
      n_valid_X <- length(W[,1])
  	}
    # fit super learner wrapper
    fit <- do.call(task$SL_wrap, list(Y = train_Y, X = train_X,
                                 newX = valid_X,
                                 obsWeights = rep(1, length(train_Y)),
                                 family = sl_control$family))

    # split up validation predictions
    if(sum(valid_idx) > 0){
      fit$pred_setA <- split(fit$pred, folds[valid_idx])
      fit$pred <- split(fit$pred, folds[valid_idx])
    }else{
      fit$pred_setA <- fit$pred
    }
    # return(fit)
    return(c(task, fit))
}

get_ate_cv_g_pred <- function(A, V, all_fit_tasks, all_fits, all_sl, folds, 
                             sl_control, return_learner_fits = TRUE,
                             learners, remove_index = NULL, compute_superlearner = TRUE){
	n <- length(A)
	train_matrix <- combn(V, V-1)
	all_out <- lapply(split(train_matrix,col(train_matrix)), function(tr){
		learner_idx <- cvma:::search_fits_for_training_folds(fits = all_fit_tasks, 
                                          y = "A", 
                                          training_folds = tr)
		if(!is.null(remove_index)){
			learner_idx <- learner_idx[-remove_index]
		}
    if(compute_superlearner){
		sl_idx <- cvma:::search_fits_for_training_folds(fits = all_sl, 
                                          y = "A", 
                                          training_folds = tr)
    }
		all_pred_setA <- Reduce(cbind,sapply(learner_idx, function(i){
			unlist(all_fits[[i]]$pred_setA)
		}, simplify = FALSE))
    if(compute_superlearner){
  		sl_pred_setA <- do.call(sl_control$ensemble_fn, args = list(pred = all_pred_setA,
  		                                                            weight = all_sl[[sl_idx]]$sl_weight))
    }else{
      sl_pred_setA <- NULL
    }
		return(list(learner_pred = all_pred_setA,
		            sl_pred = sl_pred_setA))
	})
  	
  ### !!!! NEED TO CHECK IF THIS IS RIGHT !!!! ####
  ### Looks like it might be?                  ####
  ### Assess via debug
    idx <- unlist(split(1:n, folds)[V:1], use.names = FALSE)
  	cv_learner_pred <- matrix(NA, nrow = length(A), 
  	                          ncol = ifelse(is.null(remove_index), length(learners),
  	                                        length(learners) - length(remove_index)))
  if(length(learners) > 1){
    cv_learner_pred[idx,] <- Reduce(rbind, lapply(all_out, "[[", "learner_pred"))
  }else{
    cv_learner_pred[idx,] <- Reduce(c, lapply(all_out, "[[", "learner_pred"))
  }  
  if(compute_superlearner){
	cv_sl_pred <- rep(NA, length(A))
	cv_sl_pred[idx] <- Reduce(c, lapply(all_out, "[[", "sl_pred"))
  }else{
    cv_sl_pred <- NULL
  }
    return(list(cv_learner_pred = cv_learner_pred,
                cv_sl_pred = cv_sl_pred))
}


#' estimate_nuisance
#' 
#' @importFrom future.apply future_lapply
estimate_nuisance <- function(Y, W, A, V = 5, learners, 
                              remove_learner = NULL, 
                              compute_superlearner = TRUE, 
                              which_ctmle_g = "SL.hal9002",
                      sl_control_Q = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_se",
                                   weight_fn = "weight_sl_convex",
                                   cv_risk_fn = "cv_risk_sl_r2",
                                   family = gaussian(),
                                   alpha = 0.05),
                      sl_control_g = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_nloglik",
                                   weight_fn = "weight_sl_convex",
                                   cv_risk_fn = "cv_risk_sl_r2",
                                   family = binomial(),
                                   alpha = 0.05),
                      folds){
	# get initial parameter values
    n <- length(Y)
    M <- length(learners)
    Y <- data.matrix(Y); colnames(Y) <- "Y"


    #---------------------------------
    # outcome regression
    #---------------------------------
    # all learner fitting tasks
    if(compute_superlearner){
      fold_fits <- c(V-1,V-2)      
    }else{
      fold_fits <- V-1
    }
    all_fit_tasks <- cvma:::make_fit_task_list(Ynames = "Y", learners = learners, 
                                        V = V, return_outer_sl = TRUE,
                                        fold_fits = fold_fits)

    # NOTE: could be future_lapply for parallelization
    all_fits <- future.apply::future_lapply(all_fit_tasks, FUN = get_or_fit, folds = folds, 
                              W = W, A = A, Y = Y, sl_control = sl_control_Q)

    # all super learner weight-getting tasks
    if(compute_superlearner){
      all_sl_tasks <- cvma:::make_sl_task_list(Ynames = "Y", V = V, fold_fits = c(V, V-1))
      # TO DO: I have a hunch that if future_lapply requires transferring
      #        all_fits between nodes that the communication overhead will make
      #        parallelization of this step slower than doing it sequentially 
      all_sl <- lapply(all_sl_tasks, FUN = cvma:::get_sl, 
                              Y = data.frame(Y=Y), V = V, all_fit_tasks = all_fit_tasks, 
                              all_fits = all_fits, folds = folds,
                              learners = learners, sl_control = sl_control_Q)
      all_dsl <- lapply(all_sl_tasks, FUN = cvma:::get_sl, 
                              Y = data.frame(Y=Y), V = V, all_fit_tasks = all_fit_tasks, 
                              all_fits = all_fits, folds = folds,
                              learners = learners, sl_control = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_se",
                                   weight_fn = "weight_sl_01",
                                   cv_risk_fn = "cv_risk_sl_r2",
                                   family = gaussian(),
                                   alpha = 0.05))
    }else{
      all_sl <- NULL
    }

    cv_pred <- get_ate_cv_Q_pred(Y, V, all_fit_tasks, all_fits, all_sl, folds, 
                             sl_control_Q, learners = learners, remove_index = NULL,
                             compute_superlearner = compute_superlearner)
    cv_pred_dsl <- get_ate_cv_Q_pred(Y, V, all_fit_tasks, all_fits, all_dsl, folds, 
                             sl_control_Q, learners = learners, remove_index = NULL,
                             compute_superlearner = compute_superlearner)
    # get predictions for non-CV TMLE
    pred <- Reduce(cbind, lapply(all_fits[1:length(learners)], "[[", "pred_setA"))
    if(class(pred) == "numeric"){
      pred <- matrix(pred)
    }
    if(compute_superlearner){
      sl_pred <- do.call(sl_control_Q$ensemble_fn, args = list(pred = pred, weight = 
                                                             all_sl[[1]]$sl_weight))
      dsl_pred <- do.call(sl_control_Q$ensemble_fn, args = list(pred = pred, weight = 
                                                             all_sl[[1]]$dsl_weight))
    }
    # now try without a particular learner
    if(!is.null(remove_learner)){
	    remove_index <- which(learners %in% remove_learner)    	

      all_sl_rm <- lapply(all_sl_tasks, FUN = get_rm_sl, 
                       Y = Y, V = V, all_fit_tasks = all_fit_tasks, 
                       all_fits = all_fits, folds = folds,
                       learners = learners[-remove_index], remove_index = remove_index, 
                       sl_control = sl_control_Q)
      all_dsl_rm <- lapply(all_sl_tasks, FUN = get_rm_sl, 
	                     Y = Y, V = V, all_fit_tasks = all_fit_tasks, 
	                     all_fits = all_fits, folds = folds,
	                     learners = learners[-remove_index], remove_index = remove_index, 
	                     sl_control = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_se",
                                   weight_fn = "weight_sl_01",
                                   cv_risk_fn = "cv_risk_sl_r2",
                                   family = gaussian(),
                                   alpha = 0.05))

      cv_pred_rm <- get_ate_cv_Q_pred(Y, V, all_fit_tasks, all_fits, all_sl_rm, folds, 
                               sl_control_Q, learners = learners, remove_index = remove_index,
                               compute_superlearner = compute_superlearner)	    
      cv_dpred_rm <- get_ate_cv_Q_pred(Y, V, all_fit_tasks, all_fits, all_dsl_rm, folds, 
	                             sl_control_Q, learners = learners, remove_index = remove_index,
                               compute_superlearner = compute_superlearner)

	    # predictions for sl with hal removed 
	    sl_pred_rm <- do.call(sl_control_Q$ensemble_fn, args = list(pred = pred[,-remove_index], weight = 
                                                           all_sl_rm[[1]]$sl_weight))
      dsl_pred_rm <- do.call(sl_control_Q$ensemble_fn, args = list(pred = pred[,-remove_index], weight = 
                                                           all_dsl_rm[[1]]$sl_weight))

	}

    # format everything
    Qbar_list <- list()
    if(compute_superlearner){
      # hal super learner
      Qbar_list$full_sl <- data.frame(Q0W = sl_pred[1:n],
                                     Q1W = sl_pred[(n+1):(2*n)],
                                     QAW = as.numeric(ifelse(A == 0, sl_pred[1:n], sl_pred[(n+1):(2*n)])))
      Qbar_list$full_dsl <- data.frame(Q0W = dsl_pred[1:n],
                                     Q1W = dsl_pred[(n+1):(2*n)],
                                     QAW = as.numeric(ifelse(A == 0, dsl_pred[1:n], dsl_pred[(n+1):(2*n)])))

      # cv hal super learner
      Qbar_list$cv_full_sl <- data.frame(Q0W = cv_pred$cv_sl_pred[1:n],
                                     Q1W = cv_pred$cv_sl_pred[(n+1):(2*n)],
                                     QAW = as.numeric(ifelse(A == 0, cv_pred$cv_sl_pred[1:n], cv_pred$cv_sl_pred[(n+1):(2*n)])))
      Qbar_list$cv_full_dsl <- data.frame(Q0W = cv_dpred$cv_sl_pred[1:n],
                                     Q1W = cv_dpred$cv_sl_pred[(n+1):(2*n)],
                                     QAW = as.numeric(ifelse(A == 0, cv_dpred$cv_sl_pred[1:n], cv_dpred$cv_sl_pred[(n+1):(2*n)])))
      if(!is.null(remove_learner)){
  	    # hal super learner
  	    Qbar_list$rm_sl <- data.frame(Q0W = sl_pred_rm[1:n],
                                       Q1W = sl_pred_rm[(n+1):(2*n)],
                                       QAW = as.numeric(ifelse(A == 0, sl_pred_rm[1:n], sl_pred_rm[(n+1):(2*n)])))
        Qbar_list$rm_dsl <- data.frame(Q0W = dsl_pred_rm[1:n],
  	                                   Q1W = dsl_pred_rm[(n+1):(2*n)],
  	                                   QAW = as.numeric(ifelse(A == 0, dsl_pred_rm[1:n], dsl_pred_rm[(n+1):(2*n)])))
  	    # cv hal super learner
  	    Qbar_list$cv_rm_sl <- data.frame(Q0W = cv_pred_rm$cv_sl_pred[1:n],
                                       Q1W = cv_pred_rm$cv_sl_pred[(n+1):(2*n)],
                                       QAW = as.numeric(ifelse(A == 0, cv_pred_rm$cv_sl_pred[1:n], cv_pred_rm$cv_sl_pred[(n+1):(2*n)])))
        Qbar_list$cv_rm_dsl <- data.frame(Q0W = cv_dpred_rm$cv_sl_pred[1:n],
  	                                   Q1W = cv_dpred_rm$cv_sl_pred[(n+1):(2*n)],
  	                                   QAW = as.numeric(ifelse(A == 0, cv_dpred_rm$cv_sl_pred[1:n], cv_dpred_rm$cv_sl_pred[(n+1):(2*n)])))
      }
    }
    # add other predictions
    for(i in 1:ncol(pred)){
    	Qbar_list[[learners[i]]] <- data.frame(Q0W = pred[1:n,i],
                                   Q1W = pred[(n+1):(2*n),i],
                                   QAW = as.numeric(ifelse(A == 0, pred[1:n,i], pred[(n+1):(2*n),i])))
    	Qbar_list[[paste0("cv_",learners[i])]] <- data.frame(Q0W = cv_pred$cv_learner_pred[1:n,i],
                           Q1W = cv_pred$cv_learner_pred[(n+1):(2*n),i],
                           QAW = as.numeric(ifelse(A == 0, cv_pred$cv_learner_pred[1:n,i], cv_pred$cv_learner_pred[(n+1):(2*n),i])))
    }

    #---------------------------------
    # propensity score
    #---------------------------------
	# all learner fitting tasks
    all_fit_tasks <- cvma:::make_fit_task_list(Ynames = "A", learners = learners, 
                                        V = V, return_outer_sl = TRUE,
                                        fold_fits = fold_fits)

    # NOTE: could be future_lapply for parallelization
    # NOTE: Over writes all_fits for outcome regression -- might want to change later
    all_fits <- future.apply::future_lapply(all_fit_tasks, FUN = get_ps_fit, folds = folds, 
                              W = W, A = A, sl_control = sl_control_g)

    # get CTMLE fits
    ctmle_g_fits <- get_ctmle_g_fits(all_fits = all_fits, all_fit_tasks = all_fit_tasks, 
                                 which_ctmle_g = which_ctmle_g, W = W, V = V, folds = folds)



    # all super learner weight-getting tasks
    if(compute_superlearner){
      all_sl_tasks <- cvma:::make_sl_task_list(Ynames = "A", V = V, fold_fits = c(V, V-1))
      # TO DO: I have a hunch that if future_lapply requires transferring
      #        all_fits between nodes that the communication overhead will make
      #        parallelization of this step slower than doing it sequentially 
      all_sl <- lapply(all_sl_tasks, FUN = cvma:::get_sl, 
                              Y = data.frame(A=A), V = V, all_fit_tasks = all_fit_tasks, 
                              all_fits = all_fits, folds = folds,
                              learners = learners, sl_control = sl_control_g)
      all_dsl <- lapply(all_sl_tasks, FUN = cvma:::get_sl, 
                              Y = data.frame(A=A), V = V, all_fit_tasks = all_fit_tasks, 
                              all_fits = all_fits, folds = folds,
                              learners = learners, sl_control = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_nloglik",
                                   weight_fn = "weight_sl_01",
                                   cv_risk_fn = "cv_risk_sl_r2",
                                   family = binomial(),
                                   alpha = 0.05))
    }
    cv_pred <- get_ate_cv_g_pred(A, V, all_fit_tasks, all_fits, all_sl, folds, 
                             sl_control_g, learners = learners, remove_index = NULL,
                             compute_superlearner = compute_superlearner)    
    cv_dpred <- get_ate_cv_g_pred(A, V, all_fit_tasks, all_fits, all_dsl, folds, 
                             sl_control_g, learners = learners, remove_index = NULL,
                             compute_superlearner = compute_superlearner)
    # get predictions for non-CV TMLE
    pred <- Reduce(cbind, lapply(all_fits[1:length(learners)], "[[", "pred_setA"))
    if(class(pred) == "numeric"){
      pred <- matrix(pred)
    }
    if(compute_superlearner){
      sl_pred <- do.call(sl_control_g$ensemble_fn, args = list(pred = pred, weight = 
                                                             all_sl[[1]]$sl_weight))
      dsl_pred <- do.call(sl_control_g$ensemble_fn, args = list(pred = pred, weight = 
                                                             all_dsl[[1]]$sl_weight))
    # now try without a particular learner
    if(!is.null(remove_learner)){
    remove_index <- which(learners %in% remove_learner)
	    all_sl_rm <- lapply(all_sl_tasks, FUN = get_rm_sl, 
                       Y = A, V = V, all_fit_tasks = all_fit_tasks, 
                       all_fits = all_fits, folds = folds,
                       learners = learners[-remove_index], remove_index = remove_index, 
                       sl_control = sl_control_g)
      all_dsl_rm <- lapply(all_sl_tasks, FUN = get_rm_sl, 
	                     Y = A, V = V, all_fit_tasks = all_fit_tasks, 
	                     all_fits = all_fits, folds = folds,
	                     learners = learners[-remove_index], remove_index = remove_index, 
	                     sl_control = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_nloglik",
                                   weight_fn = "weight_sl_01",
                                   cv_risk_fn = "cv_risk_sl_r2",
                                   family = binomial(),
                                   alpha = 0.05))
    	# predictions for sl with hal removed 
      sl_pred_rm <- do.call(sl_control_g$ensemble_fn, args = list(pred = pred[,-remove_index], weight = 
                                                             all_sl_rm[[1]]$sl_weight))
	    dsl_pred_rm <- do.call(sl_control_g$ensemble_fn, args = list(pred = pred[,-remove_index], weight = 
	                                                           all_dsl_rm[[1]]$sl_weight))

	   cv_pred_rm <- get_ate_cv_g_pred(A, V, all_fit_tasks, all_fits, all_sl_rm, folds, 
                               sl_control_g, learners = learners, remove_index = remove_index,
                               compute_superlearner = compute_superlearner)
     cv_dpred_rm <- get_ate_cv_g_pred(A, V, all_fit_tasks, all_fits, all_dsl_rm, folds, 
	                             sl_control_g, learners = learners, remove_index = remove_index,
                               compute_superlearner = compute_superlearner)
	   }
    }

    # oragnize predictions
    # format everything
    g_list <- list()
    if(compute_superlearner){
      # hal super learner
      g_list$full_sl <- data.frame(g1W = sl_pred, g0W = 1 - sl_pred)
      g_list$full_dsl <- data.frame(g1W = dsl_pred, g0W = 1 - dsl_pred)
      # cv hal super learner
      g_list$cv_full_sl <- data.frame(g1W = cv_pred$cv_sl_pred, g0W = 1 - cv_pred$cv_sl_pred)
      g_list$cv_full_dsl <- data.frame(g1W = cv_pred$cv_dsl_pred, g0W = 1 - cv_pred$cv_dsl_pred)
      # hal super learner
      if(!is.null(remove_learner)){
        g_list$sl <- data.frame(g1W = sl_pred_rm, g0W = 1 - sl_pred_rm)
  	    g_list$dsl <- data.frame(g1W = dsl_pred_rm, g0W = 1 - dsl_pred_rm)
  	    # cv hal super learner
        g_list$cv_sl <- data.frame(g1W = cv_pred_rm$cv_sl_pred, g0W = 1 - cv_pred_rm$cv_sl_pred)
  	    g_list$cv_dsl <- data.frame(g1W = cv_dpred_rm$cv_sl_pred, g0W = 1 - cv_dpred_rm$cv_sl_pred)
      }
    }
    # add other predictions
    for(i in 1:ncol(pred)){
    	g_list[[learners[i]]] <- data.frame(g1W = pred[ , i], g0W = 1 - pred[ , i])
    	g_list[[paste0("cv_",learners[i])]] <- data.frame(g1W = cv_pred$cv_learner_pred[ , i],
    	                                                  g0W = 1 - cv_pred$cv_learner_pred[ , i])
    }

    return(list(Qbar = Qbar_list, g = g_list, folds = folds, ctmle_g_fits = ctmle_g_fits))
}

#' Get propensity score fits in a proper format for ctmle
get_ctmle_g_fits <- function(all_fits, all_fit_tasks, which_ctmle_g, V, W, folds){
  # need to get n X K matrix of predictions from HAL when 
  # HAL is fit to the full data
  # -------------------------------------------------------
  # find full HAL fit
  full_hal_idx <- cvma:::search_fits_for_learner(training_folds = seq_len(V),
                                          fits = all_fit_tasks, 
                                          learner = which_ctmle_g,
                                          y = "A")
  g_matrix <- predict_alllambda_SL.hal9002(all_fits[[full_hal_idx]]$fit$object, newdata = W)

  # need to get n X K matrix of predictions from CV HAL 
  n <- length(W[,1])
  train_matrix <- combn(V, V-1)
  min_idx <- which(all_fits[[full_hal_idx]]$fit$object$hal_lasso$glmnet.fit$lambda == 
                    all_fits[[full_hal_idx]]$fit$object$hal_lasso$lambda.1se)
  all_out <- lapply(split(train_matrix,col(train_matrix)), function(tr){
    learner_idx <- cvma:::search_fits_for_learner(fits = all_fit_tasks, 
                                          y = "A", learner = which_ctmle_g, 
                                          training_folds = tr)
    this_g_matrix <- predict_alllambda_SL.hal9002(all_fits[[learner_idx]]$fit$object, 
                                                  newdata = W[-which(folds %in% tr),,drop=FALSE],
                                                  min_idx = min_idx)
    return(this_g_matrix)
  })

  cv_g_matrix <- matrix(NA, nrow = length(W[,1]), ncol = dim(g_matrix)[2])
  idx <- unlist(split(1:n, folds)[V:1], use.names = FALSE)
  cv_g_matrix[idx,] <- Reduce(rbind, lapply(all_out, "[[", 1))

  return(list(g_matrix = g_matrix, 
              cv_g_matrix = cv_g_matrix))
}


get_iptw <- function(Y, W, A, g, gtol = 1e-2){
  n <- length(Y)
  gAW <- rep(0, n)
  gAW[A == 1] <- g$g1W[A == 1]
  gAW[A == 0] <- g$g0W[A == 0]
  iptw_est <- mean((2*A - 1)/gAW * Y)
  return(list(est = iptw_est, se = NA))
}

get_gcomp <- function(Y, W, A, Qbar){
  return(list(est = mean(Qbar$Q1W) - mean(Qbar$Q0W), se = NA))
}

#' @param Qbar a data.frame with names QAW, Q1W, Q0W 
#' @param g a data.frame with names g1W, g0W
logistic_tmle <- function(Y, W, A, Qbar, g, gtol = 1e-2){
	l <- min(Y); u <- max(Y)
	Y_scale <- rescale(Y, l, u)
	Qbar_scale <- data.frame(apply(Qbar, 2, rescale, l=l, u=u))
	HAW <- (2*A - 1)/ifelse(A==1, g$g1W, g$g0W)
	offset <- trim_qlogis(Qbar_scale$QAW)
	fluc_dat <- data.frame(out = as.numeric(Y_scale), 
	                       off = offset, covar = as.numeric(HAW))
  suppressWarnings(
	fm <- glm(out ~ -1 + offset(off) + covar, family = binomial(),
	          data = fluc_dat)
  )
	pred_A1_dat <- data.frame(out = as.numeric(Y_scale),
	                          off = trim_qlogis(Qbar_scale$Q1W), 
	                          covar = 1/g$g1W)
	pred_A0_dat <- data.frame(out = as.numeric(Y_scale),
	                          off = trim_qlogis(Qbar_scale$Q0W), 
	                          covar = -1/g$g0W)
	Q1W_star <- re_rescale(predict(fm, newdata = pred_A1_dat, type = "response"), l, u)
	Q0W_star <- re_rescale(predict(fm, newdata = pred_A0_dat, type = "response"), l, u)
	QAW_star <- ifelse(A==1, Q1W_star, Q0W_star)
	return(data.frame(Q0W_star = Q0W_star, Q1W_star = Q1W_star,
	                  QAW_star = QAW_star))	
}

linear_tmle <- function(Y, W, A, Qbar, g, gtol = 1e-2){
	HAW_1 <- as.numeric(A==1)/g$g1W
	HAW_0 <- as.numeric(A==0)/g$g0W
	offset_1 <- Qbar$Q1W
	offset_0 <- Qbar$Q0W
	fluc_dat_1 <- data.frame(out = as.numeric(Y), 
	                       	 off = offset_1)
	fluc_dat_0 <- data.frame(out = as.numeric(Y), 
	                       	 off = offset_0)
	fm_1 <- glm(out ~ 1 + offset(off), family = gaussian(),
	          data = fluc_dat_1, weights = HAW_1)
	fm_0 <- glm(out ~ 1 + offset(off), family = gaussian(),
	          data = fluc_dat_0, weights = HAW_0)

	pred_A1_dat <- data.frame(out = as.numeric(Y),
	                          off = Qbar$Q1W)
	pred_A0_dat <- data.frame(out = as.numeric(Y),
	                          off = Qbar$Q0W)
	Q1W_star <- predict(fm_1, newdata = pred_A1_dat, type = "response")
	Q0W_star <- predict(fm_0, newdata = pred_A0_dat, type = "response")
	QAW_star <- ifelse(A==1, Q1W_star, Q0W_star)
	return(data.frame(Q0W_star = Q0W_star, Q1W_star = Q1W_star,
	                  QAW_star = QAW_star))
}

get_tmle_ate <- function(Qbar){
	return(mean(Qbar$Q1W_star) - mean(Qbar$Q0W_star))
}

get_onestep_ate <- function(Qbar, g, Y, A, W){
	HAW <- (2*A - 1)/ifelse(A==1, g$g1W, g$g0W)
	naive_ate <- mean(Qbar$Q1W - Qbar$Q0W)
	correction <- mean(HAW * (Y - Qbar$QAW))
	return(naive_ate + correction)
}

get_tmle_se <- function(Qbar, g, ate, W, A, Y){
	HAW <- (2*A - 1)/ifelse(A==1, g$g1W, g$g0W)
	Dstar <- HAW * (Y - Qbar$QAW_star) + Qbar$Q1W_star - Qbar$Q0W_star - ate
	n <- length(Y)
	return(sqrt(var(Dstar)/n))
}
get_onestep_se <- function(Qbar, g, ate, W, A, Y){
	HAW <- (2*A - 1)/ifelse(A==1, g$g1W, g$g0W)
	D <- HAW * (Y - Qbar$QAW) + Qbar$Q1W - Qbar$Q0W - ate
	n <- length(Y)
	return(sqrt(var(D)/n))
}

rescale <- function(x, l, u){
	(x-l)/(u-l)
}
re_rescale <- function(x, l, u){
	x*(u-l) + l
}
trim_qlogis <- function(x, trim = 1e-5){
	x[x < trim] <- trim; x[x > 1-trim] <- 1-trim
	qlogis(x)
}

#' Get drtmles from output of estimate nuisance
# @param Q The outcome regression formatted as output of estimate_nuisance
# @param g The propensity regression formatted as output of estimate_nuisance
get_dr_tmle <- function(W, A, Y, Q, g, folds, est_name, ...){
  Qn <- list(Q$Q0W, Q$Q1W)
  gn <- list(g$g0W, g$g1W)
  # if it's a cv estimate of nuisance, then pass in
  # the folds used to estimate
  if(grepl("cv", est_name)){
    cvFolds <- folds
  }else{
    # otherwise, don't use cv to estimate extra nuisance
    cvFolds <- 1
  }

  dr_fit <- drtmle(W = W, A = A, Y = Y, Qn = Qn, gn = gn,
                   a_0 = c(0,1), maxIter = 5, cvFolds = cvFolds, 
                   glm_Qr = "gn + I(gn^2) + I(gn^3)", 
                   glm_gr = "Qn + I(Qn^2) + I(Qn^3)",
                   verbose = FALSE)

  ci_dr_fit <- ci(dr_fit, contrast = c(-1,1))
  est <- ci_dr_fit$drtmle[1,1]
  se <- (ci_dr_fit$drtmle[1,1] - ci_dr_fit$drtmle[1,2])/qnorm(0.975)

  return(c(est = est, se = se))
}

get_dr_iptw <- function(W, A, Y, Q, g, folds, est_name, ...){
  Qn <- list(Q$Q0W, Q$Q1W)
  gn <- list(g$g0W, g$g1W)
  # if it's a cv estimate of nuisance, then pass in
  # the folds used to estimate
  if(grepl("cv", est_name)){
    cvFolds <- folds
  }else{
    # otherwise, don't use cv to estimate extra nuisance
    cvFolds <- 1
  }

  dr_fit <- adaptive_iptw(W = W, A = A, Y = Y, Qn = Qn, gn = gn,
                   a_0 = c(0,1), maxIter = 5, cvFolds = cvFolds, 
                   SL_Qr = "SL.npreg", 
                   verbose = FALSE)

  ci_dr_fit <- ci(dr_fit, contrast = c(-1,1))
  est <- ci_dr_fit$iptw_tmle[1,1]
  se <- (ci_dr_fit$iptw_tmle[1,1] - ci_dr_fit$iptw_tmle[1,2])/qnorm(0.975)

  return(c(est = est, se = se))
}

#' Matching estimators from Sekon's package
#' @param which_estimator A numeric between 1 and 4, with 1 corresponding 
#' to ps match with no bias adjust, 2 ps match with bias adjust, 3 genmatch
#' with main terms only, 4 genmatch with squared and 2-way interaction terms
#' @importFrom Matching Match GenMatch
get_matching <- function(Y, W, A, g, which_estimator, ...){
  if(which_estimator == 1){
    # standard propensity score matching 1:1 match
    rr1 <- tryCatch({
      Matching::Match(Y = Y, Tr = A, X = g$g1W, estimand = "ATE")
    }, error = function(e){
      return(list(est = NA, se = NA, message = e))
    })
  }else if(which_estimator == 2){
      # standard propensity score matching 1:1 match with bias adjustment
   rr1 <- tryCatch({ 
    Matching::Match(Y = Y, Tr = A, X = g$g1W, estimand = "ATE", BiasAdjust = TRUE)
    }, error = function(e){
      return(list(est = NA, se = NA, message = e))
    })    
  }else if(which_estimator == 3){
      # automatic matching with main terms only
    BalanceMatrix_mt <- cbind(W, g$g1W)
    gen1 <- Matching::GenMatch(Tr = A, X = W, BalanceMatrix = BalanceMatrix_mt,
                     pop.size = 1000, verbose = FALSE, print.level = 0)

    rr1 <- tryCatch({
      Matching::Match(Y = Y, Tr = A, X = BalanceMatrix_mt, Weight.matrix = gen1,
                   estimand = "ATE")
        }, error = function(e){
      return(list(est = NA, se = NA, message = e))
    })
  }else if(which_estimator == 4){
      # automatic matching with main terms, squared, terms, and interactions
    main <- model.matrix(A ~ -1 + .^2, data = W)
    pol <- model.matrix(as.formula(paste0("A ~ -1 + ", paste0("I(", colnames(W),"^2)",collapse = "+"))), data = W)
    BalanceMatrix_large <- cbind(main, pol, g$g1W)
    gen2 <- Matching::GenMatch(Tr = A, X = W, BalanceMatrix = BalanceMatrix_large,
                     pop.size = 1000, verbose = FALSE, print.level = 0)
    rr1 <- tryCatch({
          Matching::Match(Y = Y, Tr = A, X = BalanceMatrix_large, Weight.matrix = gen2,
                   estimand = "ATE")
        },  error = function(e){
      return(list(est = NA, se = NA, message = e))
  })
  }
  return(list(est = rr1$est, se = rr1$se, message = rr1$message))
}

#' @export

get_all_ates <- function(Y, W, A, V = 5, learners, 
                              remove_learner = NULL, 
                              compute_superlearner = TRUE, 
                              gtol = 1e-3, 
                      sl_control_Q = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_se",
                                   weight_fn = "weight_sl_convex",
                                   cv_risk_fn = "cv_risk_sl_r2",
                                   family = gaussian(),
                                   alpha = 0.05),
                      sl_control_g = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_nloglik",
                                   weight_fn = "weight_sl_convex",
                                   cv_risk_fn = "cv_risk_sl_r2",
                                   family = binomial(),
                                   alpha = 0.05),
                      which_dr_tmle = c("full_sl","cv_full_sl","full_dsl",
                                        "cv_full_dsl", "SL.hal9002","cv_SL.hal9002"),
                      which_dr_iptw = c("SL.hal9002","cv_SL.hal9002",
                                        "full_sl", "cv_full_sl",
                                        "full_dsl", "cv_full_dsl",
                                        "SL.gbm.caretMod"),
                      which_Match = c("SL.hal9002","SL.glm","full_sl","full_dsl"),
                      which_ctmle_g = "SL.hal9002",
                      which_ctmle_Q = c("full_sl","cv_full_sl","full_dsl","cv_full_dsl",
                                        "SL.hal9002","cv_SL.hal9002",
                                        "SL.glm")){
	# estimate nuisance
  cat("Fitting nuisance \n")
  # TO DO: make_folds function with more options, possibly from origami?
  n <- length(A)
  folds <- rep(seq_len(V), length = n)
  folds <- sample(folds)
	nuisance <- estimate_nuisance(Y = Y, W = W, A = A, V = V, learners = learners,
	                              remove_learner = remove_learner, 
	                              sl_control_Q = sl_control_Q,
	                              sl_control_g = sl_control_g,
                                compute_superlearner = compute_superlearner,
                                which_ctmle_g = which_ctmle_g, folds = folds)
  # truncate propensity estimates
  nuisance$g <- lapply(nuisance$g, function(g){
    tmp <- apply(g, 2, function(gn){ 
      gn[gn < gtol] <- gtol
      gn[gn > 1 - gtol] <- 1 - gtol
      gn
    })
    data.frame(tmp)
  })

  cat("Getting TMLEs \n")
	# get logistic tmles
	log_tmle <- mapply(Qbar = nuisance$Qbar, g = nuisance$g, logistic_tmle,
	                   MoreArgs = list(Y = Y, A = A, W = W), SIMPLIFY = FALSE)
	# get linear tmle
	lin_tmle <- mapply(Qbar = nuisance$Qbar, g = nuisance$g, linear_tmle,
	                   MoreArgs = list(Y = Y, A = A, W = W), SIMPLIFY = FALSE)
	# get tmle ate estimates
	log_tmle_ate <- lapply(log_tmle, get_tmle_ate)

	lin_tmle_ate <- lapply(lin_tmle, get_tmle_ate)
	
	onestep_ate <- mapply(Qbar = nuisance$Qbar, g = nuisance$g, get_onestep_ate,
	                   MoreArgs = list(Y = Y, A = A, W = W), SIMPLIFY = FALSE)
  # get standard errors
  log_tmle_se <- mapply(Qbar = log_tmle, g = nuisance$g, ate = log_tmle_ate,
                        get_tmle_se, MoreArgs = list(Y = Y, A = A, W = W), 
                        SIMPLIFY = FALSE) # get standard errors
  lin_tmle_se <- mapply(Qbar = lin_tmle, g = nuisance$g, ate = lin_tmle_ate,
                        get_tmle_se, MoreArgs = list(Y = Y, A = A, W = W), 
                        SIMPLIFY = FALSE)
  onestep_se <- mapply(Qbar = nuisance$Qbar, g = nuisance$g, ate = onestep_ate,
                        get_onestep_se, MoreArgs = list(Y = Y, A = A, W = W), 
                        SIMPLIFY = FALSE)

  logistic_tmle_results <- mapply(est = log_tmle_ate, se = log_tmle_se, 
                                  FUN = c, SIMPLIFY = FALSE)
  linear_tmle_results <- mapply(est = lin_tmle_ate, se = lin_tmle_se, 
                                  FUN = c, SIMPLIFY = FALSE)
  onestep_results <- mapply(est = onestep_ate, se = onestep_se, 
                                  FUN = c, SIMPLIFY = FALSE)

  cat("Getting CTMLEs \n")
  ctmle_results <- get_ctmle_results(W=W, A = A, Y = Y, V = V, folds = folds, Qs = nuisance$Qbar[which_ctmle_Q],
                                     g = nuisance$ctmle_g_fits, n = length(A), gtol = gtol)


  cat("Getting IPTW/GCOMP \n")
  # get iptw estimators
  iptw_results <- lapply(nuisance$g, get_iptw, Y = Y, A = A, W = W)

  # get gcomp estimators
  gcomp_results <- lapply(nuisance$Q, get_gcomp, Y = Y, A = A, W = W)

  # get matching estimators
  cat("Getting Matching \n")
  matching_results <- sapply(1:4, function(x){
    all_g <- lapply(nuisance$g[which_Match], get_matching, W = W, A = A, Y = Y, which_estimator = x)
    return(all_g)
  }, simplify = FALSE)

  # get tan estimator
  cat("Getting Tan & Cao \n")
  tan_results <- get_tan_est(A = A, Y = Y, W = W, 
                             family = ifelse(all(Y %in% c(0,1)), "binomial","gaussian"))
  # get cao estimator
  cao_results <- get_cao_est(A = A, Y = Y , W = W, 
                             family = ifelse(all(Y %in% c(0,1)), "binomial","gaussian"))
  cat("Getting Vermeulen \n")
  # get vermeulen estimator 1
  verm1_results <- m.biasreducedDR.identity(A = A, W = W, Y = Y)

  # get vermeulen estimators 2 for all nuisance parameters except 
  # cv versions (since no theory to back up why that would be a good idea)
  verm2_results <- lapply(nuisance$Qbar[!grepl("cv_",names(nuisance$Qbar))], 
                          data.adaptive.biasreduced.DR,
                          A = A, Y = Y, W = W)

  # get dr inference TMLEs
  cat("Getting DR-TMLEs \n")
  dr_tmle_results <- mapply(Q = nuisance$Qbar[which_dr_tmle], g = nuisance$g[which_dr_tmle], 
                    est_name = split(which_dr_tmle, 1:length(which_dr_tmle)),
                    get_dr_tmle,
                    MoreArgs = list(folds = nuisance$folds, W = W, A = A, Y = Y),
                    SIMPLIFY = FALSE)  
  
  cat("Getting DR-IPTW \n")
  dr_iptw_results <- mapply(Q = nuisance$Qbar[which_dr_iptw], g = nuisance$g[which_dr_iptw], 
                    est_name = split(which_dr_iptw, seq_along(which_dr_iptw)),
                    get_dr_iptw,
                    MoreArgs = list(folds = nuisance$folds, W = W, A = A, Y = Y),
                    SIMPLIFY = FALSE)


	return(list(logistic_tmle = logistic_tmle_results,
	            linear_tmle = linear_tmle_results,
	            onestep = onestep_results,
              dr_tmle = dr_tmle_results,
              dr_iptw = dr_iptw_results,
              ctmle = ctmle_results,
              iptw = iptw_results,
              gcomp = gcomp_results,
              match = matching_results,
              tan = tan_results,
              cao = cao_results,
              verm1 = verm1_results,
              verm2 = verm2_results))
}


get_ctmle_results <- function(Qs, g, which_ctmle_g, which_cvtmle_Q, 
                              folds, V, n, Y, A, W, gtol){
  split_folds <- split(1:n, folds)
  all_ctmles <- lapply(Qs, function(x){
    ctmle_general_fit <- ctmle::ctmleGeneral(Y = Y, A = A, W = W, Q = x[c("Q0W","Q1W")],
                                   ctmletype = 1, 
                                   gn_candidates = g$g_matrix,
                                   gn_candidates_cv = g$cv_g_matrix,
                                   folds = split_folds, V = V, gbound = gtol)
    list(est = ctmle_general_fit$est, se = sqrt(ctmle_general_fit$var.psi))
  })
  return(all_ctmles)
}


#' data.adaptive.biasreduced.DR
#' 
#' Vermeulen estimator from IJB paper
#' @export 

data.adaptive.biasreduced.DR <-
    function(A, Y, W, Qbar, 
             zeta=0.005, 
             alpha=0.05, psi.tilde=0){
        
        expit <- function(x){1/(1+exp(-x))}
        logit <- function(x){log(x/(1-x))}
        n <- length(A)
        int.cov <- cbind(rep(1,n),W)
        # colnames(dat.cov)<-
        #     paste("cov.",1:dim(cov) [2],sep="")
        # propensity score
        mler <- glm(A ~ ., data = W,family="binomial")
        ps.par1 <- fitted(mler)
        ps.par0 <- 1 - ps.par1
        # initial conditional mean outcome
        a <- min(Y) - 0.1*abs(min(Y))
        b <- max(Y) + 0.1*abs(max(Y))
        Y.star <- (Y-a)/(b-a)
        # {
        #     if(type.initQ=="par") {
        #         mley <- lm(Y.star~cov,subset=(R==1))
        #         initQ <- predict(mley,newdata=dat.cov)
        #     }
        #     else if(type.initQ=="SL"){
        #         Ym.star <- Y.star[R==1]
        #         dat.cov.m <- dat.cov[R==1,,drop=FALSE]
        #         SL.library <- c("SL.gbm.caret1","SL.step.interaction","SL.glm")
        #         initQ <- (SuperLearner(Y=Ym.star,X=dat.cov.m,
        #                               newX=dat.cov,verbose = FALSE,
        #                               SL.library=SL.library,
        #                               method="method.NNLS")$SL.predict -a)/(b-a)
        #     }else if(type.initQ=="npreg"){
        #         dat.cov.m <- dat.cov[R==1,,drop=FALSE]
        #         fm <- do.call("SL.npreg",args=list(Y=Y[R==1], X=dat.cov.m,obsWeights=rep(1,length(Y[R==1])),
        #                                          newX=dat.cov, family=gaussian()))
        #         initQ <- (fm$pred-a)/(b-a)
        #     }
        # }

        initQ1W.trunc <- ifelse(Qbar$Q1W < zeta, zeta,
                              ifelse(Qbar$Q1W > 1-zeta, 1-zeta, Qbar$Q1W))
        initQ0W.trunc <- ifelse(Qbar$Q0W < zeta, zeta,
                              ifelse(Qbar$Q0W > 1-zeta, 1-zeta, Qbar$Q0W))

        # fluctuation
        w.cov1 <- (1 - ps.par1) / ps.par1 * int.cov
        suppressWarnings(
        fluctuationQ1W <- glm(Y.star ~ -1 + ., data = w.cov1,
                           family=binomial,
                           offset=logit(initQ1W.trunc),
                           subset=(A==1))
        )
        flucQ1W <- expit(logit(initQ1W.trunc)+
                           as.vector(coef(fluctuationQ1W)%*%t(w.cov1)))                
        w.cov0 <- (1 - ps.par0) / ps.par0 * int.cov
        suppressWarnings(
        fluctuationQ0W <- glm(Y.star ~ -1 + ., data = w.cov1,
                           family=binomial,
                           offset=logit(initQ0W.trunc),
                           subset=(A==0))
        )
        flucQ0W <- expit(logit(initQ0W.trunc)+
                           as.vector(coef(fluctuationQ0W)%*%t(w.cov0)))

        # doubly robust estimator
        U <- function(R,r,Y,outcome,ps){
            outcome+as.numeric(R == r)/ps*(Y-outcome)
        }
        U1 <- U(R=A,r = 1,Y=Y.star,outcome=flucQ1W, ps=ps.par1)
        est.trunc1 <- mean(U1)
        psi1 <- (b-a)*est.trunc1+a
        U0 <- U(R=A,r =0, Y=Y.star,outcome=flucQ0W, ps=ps.par0)
        est.trunc0 <- mean(U0)
        psi0 <- (b-a)*est.trunc0+a

        # standard error
        cov_mat <- (b-a)*cov(cbind(U0,U1))
        grad <- matrix(c(-1, 1), nrow = 2)
        se.ate <- sqrt(t(grad)%*% cov_mat %*% grad / n)
        return(list(est=psi1 - psi0,
                    se=se.ate))
    }


#' m.biasreducedDR.identity
#' 
#' Code to implement Vermeulen non-data-adaptive estimator
#' @export


m.biasreducedDR.identity<-function(A,Y,W){
  n<-length(A)
  int.cov<-cbind(rep(1,n),W)
  expit<-function(x) exp(x)/(1+exp(x))
  U <- function(R,Y,X,gamma,beta){
    (R/expit(gamma%*%t(X))*(Y-beta%*%t(X))+beta%*%t(X))
  }
  min.Uint<-function(gamma){
    -mean((-R*exp(-gamma%*%t(int.cov)))+(-(1-R)*(gamma%*%t(int.cov))))
  }
  U_mat <- matrix(NA, ncol = 2, nrow = n)
  for(a in c(0,1)){
    R <- as.numeric(A == a)
    init.gamma<-coef(glm(R~.,data = W,family="binomial"))
    sol<-nlm(min.Uint, init.gamma)
    gamma.BR<-sol$estimate
    weight<-as.vector(1/exp(gamma.BR%*%t(int.cov)))
    beta.BR<-coef(lm(Y ~ -1+.,data = int.cov, subset=(R==1),weights=weight))
    U_mat[,a+1] <- U(R,Y,int.cov,gamma.BR,beta.BR)
  }
  psi1 <- mean(U_mat[,2])
  psi0 <- mean(U_mat[,1])
  cov_mat <-  cov(U_mat)
  grad <- matrix(c(-1, 1), nrow = 2)
  se.ate <- sqrt(t(grad)%*% cov_mat %*% grad / n)
  return(list(est = psi1 - psi0, se = se.ate))
}


#' getCaoEst
#' 
#' Code to implement Cao 2009 estimator
getCaoEst <- function(R,Y,cov,family){
    # fit outcome regression
    Qmod <- glm(paste0("Y ~", paste0(colnames(cov),collapse="+")),
        data = data.frame(Y, cov)[R==1,],
        family = family)

    Qn <- predict(Qmod, newdata=data.frame(Y,cov), type="response")

    # fit "enhanced propensity" regression
    negLogLik <- function(pars,R,cov){
    delta <- pars[1]; gamma <- matrix(pars[2:length(pars)],ncol=1)
    X <- data.matrix(cbind(rep(1,length(R)), cov))
    g <- 1-exp(delta + X%*%gamma)/(1 + exp(X%*%gamma))
    return(-sum(R*log(g) + (1-R)*log(1-g)))
    }

    constraint <- function(pars, R,cov){
     delta <- pars[1]; gamma <- matrix(pars[2:length(pars)],ncol=1)
     X <- data.matrix(cbind(rep(1,length(R)), cov))
     g <- 1-exp(delta + X%*%gamma)/(1 + exp(X%*%gamma))
     # first that all probs are between 0 and 1
     c1 <- as.numeric(all(g < 1) & all(g > 0))
     # now that sum of inverse weights sum to n
     c2 <- as.numeric(sum(R/g) == length(R))    
     return(c1)
    }

    suppressWarnings(
    tmp <- optim(rep(0,ncol(cov)+2), fn = negLogLik, R = R, cov = cov,
     control=list(maxit=1e2))
    )
    p <- tmp$par
    delta <- p[1]; gamma <- matrix(p[2:length(p)],ncol=1)
    X <- data.matrix(cbind(rep(1,length(R)), cov))
    gn <- 1-exp(delta + X%*%gamma)/(1 + exp(X%*%gamma))
    if(any(gn > 1) | any(gn < 0)){
        suppressWarnings(
        tmp2 <- constrOptim.nl(tmp$par, fn = negLogLik, hin = constraint, R = R, cov = cov)
        )
        p <- tmp2$par
        delta <- p[1]; gamma <- matrix(p[2:length(p)],ncol=1)
        X <- data.matrix(cbind(rep(1,length(R)), cov))
        g2 <- 1-exp(delta + X%*%gamma)/(1 + exp(X%*%gamma))
    }
    
    est <- mean(
        R*Y/gn - (R - gn)/gn * Qn
    )
    return(est)
}

#' cao.dr 
#' 
#' Compute the Cao 2009 estimator
#' @export 
get_cao_est <- function(A,Y,W, family = "gaussian"){
    # get estimate
    est0 <- getCaoEst(R=as.numeric(A == 0),Y=Y,cov=W,family=family)
    est1 <- getCaoEst(R=as.numeric(A == 1),Y=Y,cov=W,family=family)
    # get bootstrap CI
    return(list(est = est1 - est0, se = NA))
}

tan_est <- function(W, A, Y, a, family){
  Qmod <- glm(Y ~ ., data = data.frame(Y = Y, W)[A==a,], family = family)
  Qn <- predict(Qmod, type="response", newdata=W)
  gmod <- glm(as.numeric(A == a) ~ ., data=W, family="binomial")
  gn <- predict(gmod, type="response")
  
  out3 <- iWeigReg::mn.clik(y = Y, tr = A, p = gn, g = cbind(1, Qn), X = model.matrix(gmod))
  return(out3$mu)
}

#' @importFrom iWeigReg mn.clik
get_tan_est <- function(W, A, Y, family){
   psi0 <- tan_est(W, A, Y, a = 0, family = family)
   psi1 <- tan_est(W, A, Y, a = 1, family = family)
   return(list(est = psi1 - psi0, se = NA))
}




