#' @export
SL.hal9002 <- function (Y, X, newX = NULL, degrees = NULL, family = stats::gaussian(), 
    obsWeights = rep(1, length(Y)), nlambda = 200, ...) 
{
    hal_out <- fit_hal(Y = Y, X = X, degrees = degrees, yolo = FALSE,
                       family = family$family, 
                       standardize = FALSE, fit_type = "glmnet",
                       lambda =  exp(seq(3, -50, length = 200)),
                       return_lasso = TRUE)
    if (!is.null(newX)) {
        pred <- stats::predict(object = hal_out, new_data = newX)
    }
    else {
        pred <- stats::predict(object = hal_out, new_data = X)
    }
    fit <- list(object = hal_out)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- "SL.hal9001"
    return(out)
}

#' @export 
predict_alllambda_SL.hal9002 <- function(object, newdata, ...) {
  if (!is.matrix(newdata)) {
    new_data <- as.matrix(newdata)
  }

  # generate design matrix
  pred_x_basis <- hal9001:::make_design_matrix(new_data, object$basis_list)
  group <- object$copy_map[[1]]

  # reduce matrix of basis functions
  pred_x_basis <- hal9001:::apply_copy_map(pred_x_basis, object$copy_map)

  # TO DO: Probably much faster to make one matrix operation!!!
  # generate predictions
  # preds <- sapply(object$hal_lasso$lambda, function(l){
  #   predict(object$hal_lasso, newx = pred_x_basis, 
  #                  s = l)})
  n <- dim(pred_x_basis)[1]
  preds <- cbind(rep(1, n), pred_x_basis) %*% rbind(object$hal_lasso$glmnet.fit$a0, object$hal_lasso$glmnet.fit$beta)
  return(as.matrix(preds))
}



#' @export
SL.earth.cv <- function(..., nfold = 5){
  out <- SL.earth(..., nfold = nfold)
  return(out)
}

#' @export
#' 
SL.dbarts.mod = function(Y, X, newX, family, obsWeights, id,
                     sigest = NA,
                     sigdf = 3,
                     sigquant = 0.90,
                     k = 2.0,
                     power = 2.0,
                     base = 0.95,
                     binaryOffset = 0.0,
                     ntree = 200,
                     ndpost = 1000,
                     nskip = 100,
                     printevery = 100,
                     keepevery = 1,
                     keeptrainfits = T,
                     usequants = F,
                     numcut = 100,
                     printcutoffs = 0,
                     nthread = 1,
                     keepcall = T,
                     verbose = F,
                     ...) {

  SuperLearner:::.SL.require("dbarts")
  if(dim(X)[2] == 1){
    X <- X[,1]; newX <- newX[,1]
  }
  model =
    dbarts::bart(x.train = X,
                 y.train = Y,
                 # We need to pass newX in directly due to lack of prediction.
                 x.test = newX,
                 sigest = sigest,
                 sigdf = sigdf,
                 sigquant = sigquant,
                 k = k,
                 power = power,
                 base = base,
                 binaryOffset = binaryOffset,
                 weights = obsWeights,
                 ntree = ntree,
                 ndpost = ndpost,
                 nskip = nskip,
                 printevery = printevery,
                 keepevery = keepevery,
                 keeptrainfits = keeptrainfits,
                 usequants = usequants,
                 numcut = numcut,
                 printcutoffs = printcutoffs,
                 nthread = nthread,
                 keepcall = keepcall,
                 verbose = verbose)

  # TODO: there is no predict!
  #pred = predict(model, newdata = newX)
  if (family$family == "gaussian") {
    pred = model$yhat.test.mean
  } else {
    # No mean is provided for binary Y :/
    pred = colMeans(pnorm(model$yhat.test))
  }

  fit = list(object = model)
  class(fit) = c("SL.dbarts")
  out = list(pred = pred, fit = fit)
  return(out)
}

#' SL.svmLinear.caretMod
#' Uses SL.caretMod to train a linear svm with 15 choices of 
#' tuning parameter
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export

SL.svmLinear.caretMod <- function(...,method="svmLinear",tuneLength = 15, 
                                  trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}

#' SL.rpart.caretMod 
#' 
#' Uses SL.caretMod with five-fold CV to train regression tree 
#' with ten choices of tuning parameter.
#' 
#' @param method Character describing what algorithm to train
#' @param tuneLength The number of tuning parameter combinations
#' @param trControl Passed to \code{caret::train}
#' @param ... Other arguments passed to SL.caretMod
#' 
#' @export

SL.rpart.caretMod <- function(...,method="rpart",tuneLength = 10, trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength, trControl = trControl)
}

#' SL.rf.caretMod
#' 
#' Uses SL.caretMod with five-fold CV to train random forest
#' with ten choices of tuning parameter.
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export

SL.rf.caretMod <- function (..., method = "rf", tuneLength = 10, 
                            trControl = caret::trainControl(method = "cv", number = 5)) 
{
    SL.caretMod(..., method = method, tuneLength = tuneLength, trControl = trControl)
}

#' SL.randomGLM.caretMod
#' Uses SL.caretMod with five-fold CV to train a random GLM
#' with ten choices of tuning parameter.
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export 
SL.randomGLM.caretMod <- function(...,method="randomGLM",tuneLength=10, 
                                  trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}


#' SL.nnet.caretMod
#' 
#' Uses SL.caretMod to train with five-fold CV to train a neural network
#' with twenty choices of tuning parameter.
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export 
SL.nnet.caretMod <- function(...,method="nnet", tuneLength=20, trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}

#' SL.glmnet.caretMod
#' 
#' Uses SL.caretMod with five-fold CV to train a glmnet regression
#' with ten choices of tuning parameter.
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export 
SL.glmnet.caretMod <- function(...,method="glmnet", tuneLength=10, trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}


#' SL.gbm.caretMod
#' 
#' Uses SL.caretMod to train a gbm with 10 choices of tuning parameters and 5-fold CV
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export

SL.gbm.caretMod <- function (..., method = "gbm", tuneLength = 5, 
                           trControl = caret::trainControl(method = "cv", number = 5)) 
{
    SL.caretMod(..., method = method, tuneLength = tuneLength, trControl = trControl)
}

#' SL.gamSpline.caretMod
#' 
#' Uses SL.caretMod with 10 choices of tuning parameters and 5-fold CV
#' to train a generalized additive model.
#' 
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#'  
#' @export 
SL.gamSpline.caretMod <- function(...,method="gamSpline",tuneLength=10, 
                                  trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}

#' SL.caretMod 
#' 
#' A modification of the SL.caret function from the SuperLearner package that
#' suppresses some output when method = "gbm" or "nnet". See \code{?SL.caret}
#' for more information.
#' 
#' @param Y Training outcomes
#' @param X Training predictors 
#' @param newX Test set predictors 
#' @param obsWeights Weights for the observations
#' @param method Character describing what algorithm to train
#' @param tuneLength The number of tuning parameter combinations
#' @param trControl Passed to \code{caret::train}
#' @param family Character indicating family argument (ignored)
#' @param ... Other arguments (not currently used)
#' 
#' @importFrom caret train
#' @importFrom stats predict
#' 
#' @export

SL.caretMod <- function (Y, X, newX, family, obsWeights, method, tuneLength = 10, 
                       trControl = caret::trainControl(method = "cv", number = 5, verboseIter = FALSE), 
                       ...) 
{
  if (length(unique(Y))>2){
    if(is.matrix(Y)) Y <- as.numeric(Y)
    metric <- "RMSE"
    if(method=="gbm"){
      suppressWarnings(
        # pass verbose==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,verbose=FALSE)
      )
    }else if(method=="nnet"){
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,trace=FALSE)
      )
    }else{
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl)
      )
    }
    pred <- stats::predict(fit.train, newdata = newX, type = "raw")
  }
  if (length(unique(Y))<=2) {
    metric <- "Accuracy"
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    if(method=="gbm"){
      suppressWarnings(
        # pass verbose==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl, verbose = FALSE)
      )
    }else if(method=="nnet"){
      suppressWarnings(
        # pass trace==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,trace=FALSE)
      )
    }else{
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl)
      )
    }
    pred <- stats::predict(fit.train, newdata = newX, type = "prob")[,2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret")
  return(out)
}



