
#' Une belle documentation pour RSF_regression
#'
#' ma description
#'
#' @param y_var A character string which gives the name of the right-censored variable
#' @param delta_var A character string which gives the name of the event variable (this variable should be binary : 1 = non censored, 0 = censored)
#' @param x_vars A vector of character strings which gives the names of the explanatory variables
#' @param data_train A data.frame of training observations. It must contain the column names \code{y_var},
#' \code{delta_var} and \code{x_vars}
#' @param data_test A data.frame of testing observations (default = \code{NULL})
#' @param phi A function to be applied to \code{y_var}
#' @param phi.args A list of additional parameters for the function \code{phi} (default = NULL). See \emph{Examples} for a use case
#' @param max_time A real number giving a threshold for \code{y_var} (default = \code{NULL}). If \code{NULL}, then max_time is
#' set to the maximum non censored observation of \code{y_var} among the training set
#' @param RSF_object A boolean which indicates if the Cox model fitted to the training data should be returned (default = \code{TRUE})
#' @param eval_methods A vector of character strings which gives the methods that should be used for the evaluation of the
#' model (default = \code{c("concordance","weighted")}).
#' Possible choices are "concordance", "weighted", "group" and "single". Multiple choices are possible.
#' See \emph{Details - Evaluation criteria} for more information
#' @param v_bandwidth A vector of real numbers for the bandwidths to use for the model evaluation if \code{"group"} is used
#' as an \code{eval_method} (default = \code{c(20)}). Only used if \code{"group"} is used as \code{eval_method}.
#' See \emph{Details - Evaluation criteria} for more information
#' @param types_weights_eval A vector of character strings which gives the types of weights to be used for IPCW
#' in the model evaluation (default = \code{c("KM")} (Kaplan Meier)).
#' Possible choices are "KM", "Cox", "RSF" and "0_1". See \emph{Details - Evaluation criteria} for more information
#' @param max_ratio_weights_eval A real number which gives the maximum admissible ratio for the IPC weights (default = 1000).
#' See \emph{Details - Evaluation criteria} for more information
#' @param mat_weights A matrix to provide handmade IPC weights for the model evaluation (default = \code{NULL}).
#' \code{mat_weights} should satisfied \code{nrow(mat_weights) = nrow(data_train) + nrow(data_test)} and should give the
#' in columns (multiple columns are possible). Column names of \code{mat_weights} may be used to specify names
#' for the provided weights (by default names will be "w1), "w2", ...
#' @param y_non_censored_var A character string which gives the name of the non censored \code{y_var} (default = NULL).
#' To be used only in the context os simulated data.
#' @param ... Additional parameter that may be pass to the \code{\link[survival]{coxph}} function (package \emph{survival})
#'
#'
#' @details
#' \itemize{ %je peux utiliser enumerate si je souhaite mettre des numero
#' \item \emph{Evaluation criteria}
#'
#' balbla
#' \item rezre erae}
#'
#'
#'
#'
#' @return A list with the following elements :
#' \item{predicted}{A vector of the predicted values for phi(T)}
#' \item{survival}{A dataframe (matrix) of the predicted survival curves given by the Cox model}
#' \item{time_points}{A vector of the time points where the survival curves are evaluated}
RSF_regression = function(y_var,
                          delta_var,
                          x_vars,
                          data_train,
                          data_test = NULL,
                          phi = function(x){x},
                          phi.args = list(),
                          max_time = NULL,
                          RSF_object = T,
                          eval_methods = c("concordance","weighted"),
                          v_bandwidth = 20,
                          types_weights_eval = c("KM"),
                          max_ratio_weights_eval = 20,
                          mat_weights = NULL,
                          y_non_censored_var = NULL,
                          ...){

  # Preprocessing of the arguments & data

  # preprocessing() # this doesn't work for the moment

  # column names of mat_weights should be explicit
  if(!is.null(mat_weights) & is.null(colnames(mat_weights))) colnames(mat_weights) = paste0("w",1:ncol(mat_weights))

  eval_methods <- match.arg(as.character(eval_methods), c("concordance","single", "group", "weighted"), several.ok = T)
  types_weights_eval = match.arg(as.character(types_weights_eval), c("KM", "Cox", "RSF", "0_1"), several.ok = T)

  if (is.null(y_non_censored_var)) {
    phi_non_censored_name = NULL
  } else {
    phi_non_censored_name = "phi_non_censored"
  }

  if (!is.null(data_test)){
    data = rbind(data_train[,c(y_var, delta_var, x_vars, y_non_censored_var)],
                 data_test[,c(y_var, delta_var, x_vars, y_non_censored_var)])
    data$is_train = c(rep(1, nrow(data_train)), rep(0, nrow(data_test)))
  } else {
    data = data_train[,c(y_var, delta_var, x_vars, y_non_censored_var)]
    data$is_train = 1
  }

  if (is.null(max_time)){max_time = max(data_train[which(data_train[, delta_var] == 1), y_var])}

  data$y_prime = pmin(data[,y_var], max_time)
  data$delta_prime = 1 * ((data[,delta_var] != 0) | (data[,y_var] >= max_time))
  data$phi = sapply(X = 1:length(data$y_prime),
                    FUN = function(i){do.call(phi, c(list(x=data$y_prime[i]), phi.args))})
  if(!is.null(y_non_censored_var)){
    data$phi_non_censored = sapply(X = 1:nrow(data),
                                   FUN = function(i){do.call(phi, c(list(x=pmin(data[,y_non_censored_var], max_time)[i]), phi.args))})
  }


  # Computation of the weitghts if not provided
  if (is.null(mat_weights)){
    mat_weights_train = matrix(rep(0, length(types_weights_eval) * sum(data$is_train == 1) ), ncol = length(types_weights_eval))
    colnames(mat_weights_train) = types_weights_eval
    if (!is.null(data_test)){
      mat_weights_test = matrix(rep(0, length(types_weights_eval) * sum(data$is_train == 0) ), ncol = length(types_weights_eval))
      colnames(mat_weights_test) = types_weights_eval
    }
    for (j in 1:length(types_weights_eval)){
      mat_weights_train[,j] = make_weights(data = data[data$is_train == 1, ],
                                           y_name = "y_prime",
                                           delta_name = "delta_prime",
                                           y_name2 = y_var,
                                           delta_name2 = delta_var,
                                           type = types_weights_eval[j],
                                           max_ratio_weights = 1000,
                                           x_vars = x_vars,
                                           censoring_model_object = FALSE)$weights
      if (!is.null(data_test)){
        mat_weights_test[,j] = make_weights(data = data[data$is_train == 0, ],
                                            y_name = "y_prime",
                                            delta_name = "delta_prime",
                                            y_name2 = y_var,
                                            delta_name2 = delta_var,
                                            type = types_weights_eval[j],
                                            max_ratio_weights = 1000,
                                            x_vars = x_vars,
                                            censoring_model_object = FALSE)$weights
      }
    }
  }

  # Build mat_weights_train & mat_weights_test if mat_weights provided
  if (!is.null(mat_weights)){
    mat_weights_train = mat_weights[1:nrow(data_train),]
    if (!is.null(data_test)){
      mat_weights_test = mat_weights[(nrow(data_train)+1):(nrow(data)),]
    }
  }

  # Thresholding of the weights_eval
  ## train
  mat_weights_train = apply(X = mat_weights_train, MARGIN = 2,
                            FUN = function(x){
                              x = pmin(x, min(x[x > 0]) * max_ratio_weights_eval)
                              x = x / sum(x)
                            })
  ## test
  if (!is.null(data_test)){
    mat_weights_test = apply(X = mat_weights_test, MARGIN = 2,
                             FUN = function(x){
                               x = pmin(x, min(x[x > 0]) * max_ratio_weights_eval)
                               x = x / sum(x)
                             })
  }

  # Build train & test
  data_train = data[data$is_train == 1,]
  if (!is.null(data_test)){
    data_test = data[data$is_train == 0,]
  }


  # Calibration of RSF model
  formula = stats::as.formula(paste0("Surv(", y_var, ",", delta_var," ) ~ ."))
  ntime = seq(from = 0, to = max_time * 1.05, length.out = 100)

  rfSRC = randomForestSRC::rfsrc(formula = formula ,
                                 data = data_train[,c(y_var, delta_var, x_vars)],
                                 forest = T,
                                 ntime = ntime,
                                 ...)

  overfitted_predictions_direct_RSF =
    (cbind(1, rfSRC$survival[,which(rfSRC$time.interest < max_time)]) -
       cbind(rfSRC$survival[,which(rfSRC$time.interest < max_time)],0)) %*%
    sapply(X = 1:length(c(rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)),
           FUN = function(i){do.call(phi,
                                     c(list(x=c(rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)[i]),
                                       phi.args))})

  # Performances on train test
  list_criteria_train = eval_model(predictions = overfitted_predictions_direct_RSF,
                                   data = data_train,
                                   phi_name = "phi",
                                   y_name = "y_prime",
                                   delta_name = "delta_prime",
                                   max_time = max_time,
                                   eval_methods = eval_methods,
                                   phi = phi,
                                   phi.args = phi.args,
                                   mat_weights = mat_weights_train,
                                   phi_non_censored_name = phi_non_censored_name,
                                   v_bandwidth = v_bandwidth)

  if(!is.null(data_test)){

    # Predictions on test set
    test_predictions_surv_curves_direct_RSF =
      randomForestSRC::predict.rfsrc(rfSRC, data_test[,x_vars])$survival

    test_predictions_direct_RSF =
      ( cbind(1,test_predictions_surv_curves_direct_RSF[,which(rfSRC$time.interest < max_time)]) -
          cbind(test_predictions_surv_curves_direct_RSF[,which(rfSRC$time.interest < max_time)],0) ) %*%
      sapply(X = 1:length(c(rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)),
             FUN = function(i){do.call(phi,
                                       c(list(x=c(rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)[i]),
                                         phi.args))})

    # Performances on test set
    list_criteria_test = eval_model(predictions = test_predictions_direct_RSF,
                                    data = data_test,
                                    phi_name = "phi",
                                    y_name = "y_prime",
                                    delta_name = "delta_prime",
                                    max_time = max_time,
                                    eval_methods = eval_methods,
                                    phi = phi,
                                    phi.args = phi.args,
                                    mat_weights = mat_weights_test,
                                    phi_non_censored_name = phi_non_censored_name,
                                    v_bandwidth = v_bandwidth)
  }

  result = list(
    predicted_train = as.vector(overfitted_predictions_direct_RSF),
    list_criteria_train = list_criteria_train,
    data_train = data_train[,c(y_var, delta_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)],
    mat_weights_train = mat_weights_train,
    max_time = max_time,
    phi = phi,
    phi.args = phi.args,
    x_vars = x_vars,
    censoring_rate_with_threshold = sum(data$delta_prime == 0) / nrow(data)
  )
  if (!is.null(data_test)){
    result$predicted_test = as.vector(test_predictions_direct_RSF)
    result$list_criteria_test = list_criteria_test
    result$data_test = data_test[,c(y_var, delta_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)]
    result$mat_weights_test = mat_weights_test
  }
  if (RSF_object){
    result$RSF_object = rfSRC
    result$survival_train = cbind(1,rfSRC$survival[,which(rfSRC$time.interest < max_time)], 0)
    result$time_points = c(0,rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)
    if (!is.null(data_test)){
      result$survival_test =
        cbind(1,test_predictions_surv_curves_direct_RSF[,which(rfSRC$time.interest < max_time)],0)
    }
  }
  return(result)
}

predict_RSF_regression = function(object, newdata){
  if (is.null(object$RSF_object)) stop("to use predict_RSF_regression on a RSF_regression object,
                                       you shoud specify RSF_object = TRUE in the call of RSF_regression")

  predictions_surv_curves =
    randomForestSRC::predict.rfsrc(object$RSF_object,
                                   newdata[,object$x_vars])$survival

  # compute predicted values for phi
  predictions =
    ( cbind(1,predictions_surv_curves[,which(object[["RSF_object"]][["time.interest"]] < object$max_time)]) -
        cbind(predictions_surv_curves[,which(object[["RSF_object"]][["time.interest"]] < object$max_time)],0) ) %*%
    sapply(X = 1:(sum(object[["RSF_object"]][["time.interest"]] < object$max_time) + 1),
           FUN = function(i){do.call(object$phi,
                                     c(list(x=c(object[["RSF_object"]][["time.interest"]][which(object[["RSF_object"]][["time.interest"]] < object$max_time)],
                                                object$max_time)[i]),
                                       object$phi.args))})

  return(list(predicted = as.vector(predictions),
              survival = cbind(1,predictions_surv_curves[,which(object[["RSF_object"]][["time.interest"]] < object$max_time)],0),
              time_points = c(0,object[["RSF_object"]][["time.interest"]][which(object[["RSF_object"]][["time.interest"]] < object$max_time)], object$max_time)
  ))
}


