
#' @title Fit a forest of Relative Risk Tree (RRT) and compute predictions of expectations
#' for \code{phi}\eqn{(T)}
#'
#' @description \code{rrt_reg} is a benchmark model we use in [Gerber et al. (2018)].
#' To model the variable \code{phi}\eqn{(T)}, where \eqn{T} is a right censored time and
#' \code{phi} is a given function, we first
#' fit a forest RTT model to the data to estimate the survival function
#' of \eqn{T} given the covariates. Then,
#' we deduce an estimator of \code{phi}\eqn{(T)} by integration of the function \code{phi} with respect
#' to the estimated survival function. Different methods are available to assess the quality of fit of
#' \code{rrt_reg}. \code{rrt_reg} is a wrapper for the \code{\link[rpart]{rpart}}
#' function\cr \cr
#' The notations we use are :
#' \itemize{
#' \item \eqn{C} : Censoring variable
#' \item \eqn{Y = min(T, C)}
#' \item \eqn{\delta = 1_{T \le C}}  (delta)
#' }
#'
#'
#' @param y_var A character string which gives the name of the \eqn{Y} variable
#'
#' @param delta_var A character string which gives the name of the \eqn{\delta} variable (this variable
#' should be binary : 1 = non censored, 0 = censored)
#'
#' @param x_vars A vector of character strings which gives the names of the explanatory variables
#'
#' @param train A data.frame of training observations. It must contain the column names \code{y_var},
#' \code{delta_var} and \code{x_vars}
#'
#' @param test A data.frame of testing observations (default = \code{NULL})
#'
#' @param phi A function to be applied to \code{y_var}
#'
#' @param phi.args A list of additional parameters for the function \code{phi} (default = NULL).
#' See \emph{Examples} for a use case
#'
#' @param max_time A real number giving a threshold for \code{y_var} (default = \code{NULL}).
#' If \code{NULL}, then \code{max_time} is set to the maximum non censored observation
#' of \code{y_var} among the training set. We note \eqn{T' = min(T,} \code{max_time}\eqn{)}
#'
#' @param rrt_obj A boolean which indicates if the relative risk forest object fitted to the
#' training data should
#' be returned (default = \code{TRUE})
#'
#' @param ev_methods A vector of character strings which gives the methods that should be
#' used for the evaluation of the model (default = \code{c("concordance","weighted")}).
#' Possible choices are \code{"concordance"}, \code{"weighted"} and \code{"group"}. Multiple choices are possible.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param bandwidths A vector of real numbers for the bandwidths to use for the model
#' evaluation if \code{"group"} is used
#' as an \code{ev_methods} (default = \code{NULL} : set to 50 if \code{"group"} in \code{ev_methods}).
#' Only used if \code{"group"} is used as \code{ev_methods}.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param types_w_ev A vector of character strings which gives the types of weights to be used for IPCW
#' (Inverse Probability of Censoring Weighting) in the model evaluation (default = \code{c("KM")} (Kaplan Meier)).
#' Possible choices are \code{"KM"}, \code{"Cox"}, \code{"RSF"} and \code{"unif"}. Set to \code{colnames(mat_w)}
#' if \code{mat_w} is provided. See
#' \emph{Details - Evaluation criteria} for more information
#'
#' @param max_w_ev A real number which gives the maximum admissible ratio
#' for the IPC weights (default = 1000).
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param mat_w A matrix to provide handmade IPC weights for the model evaluation (default = \code{NULL}).
#' \code{mat_w} should satisfied \code{nrow(mat_w) = nrow(train) + nrow(test)} and a column
#' should correspond to a type of weights (multiple columns are possible).
#' Column names of \code{mat_w} may be used to specify names
#' for the provided weights (by default names will be "w1", "w2", ...)
#'
#' @param ntree A positive integer which gives the number of trees to grow
#' in the forest (default = \code{100}).
#'
#' @param minleaf A positive integer indicating the minimum number of observations
#' that must be present in a (terminal) leaf (default = \code{5}).
#'
#' @param maxdepth A positive integer indicating the maximum number of layers
#' in individual trees (default = \code{6}).
#'
#' @param mtry A positive integer indicating the number of random variables
#' to draw at each node for the split selection (default = \code{NULL}).
#' If \code{NULL}, \code{mtry} is set to \code{floor(sqrt(length(x_vars)))} by default.
#'
#' @param y_no_cens_var A character string which gives the name of the non censored \code{y_var} (default = NULL).
#' To be used only in the context of simulated data where full about is available.
#'
#' @param ... Additional parameter that may be pass to the \code{\link[rpart]{rpart}}
#' function (package \emph{rpart})
#'
#'
#' @details
#' \itemize{
#' \item \emph{Evaluation criteria}
#'
#' The quality of fit may be assess throught three different criteria. There are two main
#' criteria : \code{"weighted"} and \code{"concordance"}, and one criteria that is expimental : \code{"group"}.
#'
#' \itemize{
#'
#' \item \code{"weighted"} : the weighed criteria is described in §? of [Gerb. et al.]. This criteria aims
#' to estimate the quadratic error of the model in the context of right censoring. It has the form
#' \eqn{\sum_i W_i  (y_i - \hat{y}_i)^2} where \eqn{(y_i)i} are the censored targets of the model, \eqn{(W_i)i}
#' are the IPC weights, and \eqn{(\hat{y}_i)i} are the predictions made.
#'
#' The \code{types_w_ev} argument allows the use of four kinds of IPC weights :
#' \code{"KM"}, \code{"Cox"}, \code{"RSF"} and \code{"unif"}. The first three types of weights correspond
#' to different ways to estimate the survival function of the censoring. On the other hand, \code{"unif"}
#' corresponds to \eqn{W_i = 1} for all i.
#'
#' Since the IPC weights may take too large values in some situation, \code{max_w_ev} allows
#' to threshold the weights \eqn{(W_i)i} so that the ratio between the largest and the smallest weights doesn't
#' exceed \code{max_w_ev}. If \eqn{W_max} is the maximum weight, considered weights are then
#' \eqn{min(W_i, W_max) / ( \sum_i min(W_i, W_max) ) }
#'
#' You can also manually provide weights to be used for IPCW with the argument \code{mat_w} (those
#' weights will also be threshold w.r.t. \code{max_ration_weights_eval}). \code{mat_w} should be a
#' matrix satisfying \code{nrow(mat_w) = nrow(train) + nrow(test)}, any numbers of
#' columns may be provided and then each column of the matrix corresponds to a type of weights.
#' Columns name of the matrix are taken as names for the different types of weights. If there is no
#' column name, then default names are "w1", "w2', ...
#'
#'
#' \item \code{"concordance"} : The concordance is a classical measure of performance when modelling
#' censored variables. It indicates if the order of the predicted values of the model is similar to
#' the order of the observed values. The concordance generalizes the Kendall tau to the censored case.
#'
#'
#' \item \code{"group"} : This is an experimental criteria that we didn't mention in [Gerb. et al.]. Here,
#' the idea to take the censoring into account is to measure errors given groups of observations and
#' not single observations. First, the test sample is ordered w.r.t. the predicted values of the
#' model. Second, respecting this order the test sample is splited into groups of size \code{v_bandwidht[1]}, or
#' \code{v_bandwidht[2]}, etc ... (each bandwidth in \code{v_bandwidht} corresponds to a different
#' score \code{"group"}).
#'
#' Then, inside a group, an estimator of the survival function of \eqn{T} may be obtained by
#' Kaplan Meier, and we can deduce an estimator of \code{phi}\eqn{(T)} by integration. This estimator
#' of \code{phi}\eqn{(T)} may be
#' viewed as an "empirical" value of \code{phi}\eqn{(T)} inside the group.
#'
#' On the other other hand, each observation of a group is associated to a prediction of \code{phi}\eqn{(T)}.
#' The prediction for the group may be defined as the mean of the predictions of the observations
#' inside the group.
#'
#' The resulting group scores correspond to the classical quality of fit criteria (e.g. quadratic error,
#' Kendall tau, Gini index) but taken on groups and not on single observations.
#'
#' The \code{"group"} criteria is interesting in the context of big database, in which sufficient
#' number of groups are available. This criteria has a high variance if applied to small test sample.
#'
#' }
#' }
#'
#' @return A list with the following elements :
#' \item{pred_train}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' (with \eqn{T' = min(T, } \code{max_time}\eqn{)}) for the observations of the train set}
#' \item{pred_test}{The vector of the predicted values for
#' \code{phi}\eqn{(T')} for the observations of the test set (require \code{test} != \code{NULL})}
#' \item{perf_train}{The list with the values for the evaluation criteria computed on the train
#' set}
#' \item{perf_test}{The list with the values for the evaluation criteria computed on the test
#' set (require \code{test} != \code{NULL}))}
#' \item{surv_train}{The matrix which contains the estimated values of the survival curves at
#' \code{time_points} (with the Cox model), for the observations of the train set}
#' \item{surv_test}{The matrix which contains the estimated values of the survival curves at
#' \code{time_points}, for the observations of the test set (require \code{test} != \code{NULL})}
#' \item{time_points}{The vector of the time points where the survival curves are evaluated}
#'
#' \item{mat_w_train}{The matrix which contains the values of the weights used for the
#' \code{"weighted"} criteria, for the observations of the train set}
#'
#' \item{mat_w_test}{The matrix which contains the values of the weights used for the
#' \code{"weighted"} criteria, for the observations of the test set}
#
#' \item{n_w_ev_modif_train}{The vector giving the number of train weights modified due to
#' \code{max_w_ev}}
#'
#' \item{n_w_ev_modif_test}{The vector giving the number of test weights modified due to
#' \code{max_w_ev}}
#'
#' \item{rrt_obj}{The relative risk forest object}
#'
#' \item{max_time}{The real number giving the threshold used by the model}
#'
#' \item{cens_rate}{The real number giving the rate of censoring
#' of \eqn{T'}, computed on the concatenation of \code{train} and \code{test}}
#'
#' \item{train}{The data.frame of the train data provided as arguments, plus columns :
#' \eqn{Y' = min(Y,} \code{max_time} \eqn{)}, \eqn{\delta' = 1_{T' \le C}}
#' and \code{phi}\eqn{(T')} }
#' \item{test}{The data.frame of the test data provided as arguments, plus columns :
#' \eqn{Y' = min(Y,} \code{max_time} \eqn{)}, \eqn{\delta' = 1_{T' \le C}}
#' and \code{phi}\eqn{(T')} }
#'
#' \item{phi}{See \emph{Argument}}
#'
#' \item{phi.args}{See \emph{Argument}}
#'
#' \item{x_vars}{See \emph{Argument}}
#'
#' \item{max_w_ev}{See \emph{Argument}}
#'
#'
#' @references Gerber, G., Le Faou, Y., Lopez, O., & Trupin, M. (2018). \emph{The impact of churn on
#' prospect value in health insurance, evaluation using a random forest under random censoring.}
#' \url{https://hal.archives-ouvertes.fr/hal-01807623/}
#'
#' @seealso \code{\link[rpart]{rpart}}, \code{\link{predict_rrt_reg}}
#'
#' @export
#'
#'
rrt_reg = function(y_var,
                   delta_var,
                   x_vars,
                   train,
                   test = NULL,
                   phi = function(x){x},
                   phi.args = list(),
                   max_time = NULL,
                   rrt_obj = T,
                   ev_methods = c("concordance","weighted"),
                   bandwidths = NULL,
                   types_w_ev = c("KM"),
                   max_w_ev = 20,
                   mat_w = NULL,
                   y_no_cens_var = NULL,

                   # param. for RF
                   ntree = 100,
                   minleaf = 5,
                   maxdepth = 6,
                   mtry = NULL,
                   ...){

  # Preprocessing of the arguments & data

  # column names of mat_w should be explicit
  if(!is.null(mat_w) & is.null(colnames(mat_w))) colnames(mat_w) = paste0("w",1:ncol(mat_w))

  ev_methods <- match.arg(as.character(ev_methods), c("concordance", "group", "weighted"), several.ok = T)
  if (is.null(bandwidths) & ("group" %in% ev_methods)) bandwidths = 50

  if (is.null(mat_w)){
    types_w_ev = match.arg(as.character(types_w_ev), c("KM", "Cox", "RSF", "unif"), several.ok = T)
  }

  if(is.null(mtry)){mtry = floor(sqrt(length(x_vars)))}

  if (is.null(y_no_cens_var)) {
    phi_non_censored_name = NULL
  } else {
    phi_non_censored_name = "phi_non_censored"
  }

  if (!is.null(test)){
    data = rbind(train[,c(y_var, delta_var, x_vars, y_no_cens_var)],
                 test[,c(y_var, delta_var, x_vars, y_no_cens_var)])
    data$is_train = c(rep(1, nrow(train)), rep(0, nrow(test)))
  } else {
    data = train[,c(y_var, delta_var, x_vars, y_no_cens_var)]
    data$is_train = 1
  }

  if (("group" %in% ev_methods) & (nrow(data) < 500)){
    bandwidths = pmin(bandwidths, ifelse(!is.null(test), nrow(test), nrow(train)))
    stop("group performance criteria must not be accurate because it
         needs more observations to converge")
  }
  if (is.null(max_time)){max_time = max(train[which(train[, delta_var] == 1), y_var])}

  data$y_prime = pmin(data[,y_var], max_time)
  data$delta_prime = 1 * ((data[,delta_var] != 0) | (data[,y_var] >= max_time))
  data$phi = sapply(X = 1:length(data$y_prime),
                    FUN = function(i){do.call(phi, c(list(x=data$y_prime[i]), phi.args))})
  if(!is.null(y_no_cens_var)){
    data$phi_non_censored = sapply(X = 1:nrow(data),
                                   FUN = function(i){do.call(phi, c(list(x=pmin(data[,y_no_cens_var], max_time)[i]), phi.args))})
  }


  # Computation of the weitghts if not provided
  if (is.null(mat_w)){
    mat_w_train = matrix(rep(0, length(types_w_ev) * sum(data$is_train == 1) ), ncol = length(types_w_ev))
    colnames(mat_w_train) = types_w_ev
    if (!is.null(test)){
      mat_w_test = matrix(rep(0, length(types_w_ev) * sum(data$is_train == 0) ), ncol = length(types_w_ev))
      colnames(mat_w_test) = types_w_ev
    }
    for (j in 1:length(types_w_ev)){
      mat_w_train[,j] = make_weights(data = data[data$is_train == 1, ],
                                     y_name = "y_prime",
                                     delta_name = "delta_prime",
                                     y_name2 = y_var,
                                     delta_name2 = delta_var,
                                     type = types_w_ev[j],
                                     max_ratio_weights = 1000,
                                     x_vars = x_vars,
                                     cens_mod_obj = FALSE)$weights
      if (!is.null(test)){
        mat_w_test[,j] = make_weights(data = data[data$is_train == 0, ],
                                      y_name = "y_prime",
                                      delta_name = "delta_prime",
                                      y_name2 = y_var,
                                      delta_name2 = delta_var,
                                      type = types_w_ev[j],
                                      max_ratio_weights = 1000,
                                      x_vars = x_vars,
                                      cens_mod_obj = FALSE)$weights
      }
    }
  }

  # Build mat_w_train & mat_w_test if mat_w provided
  if (!is.null(mat_w)){
    if (identical(types_w_ev, "KM")){
      types_w_ev = colnames(mat_w)
    }
    mat_w_train = as.matrix(mat_w[1:nrow(train), types_w_ev])
    if (!is.null(test)){
      mat_w_test = as.matrix(mat_w[(nrow(train)+1):(nrow(data)), types_w_ev])
    }
  }

  # Thresholding of the weights_eval
  ## train

  n_w_ev_modif_train = apply(X = mat_w_train, MARGIN = 2,
                             FUN = function(x){
                               x = sum(x > min(x[x > 0]) * max_w_ev)
                             })

  mat_w_train = apply(X = mat_w_train, MARGIN = 2,
                      FUN = function(x){
                        x = pmin(x, min(x[x > 0]) * max_w_ev)
                        x = x / sum(x)
                      })
  ## test
  if (!is.null(test)){

    n_w_ev_modif_test = apply(X = mat_w_test, MARGIN = 2,
                              FUN = function(x){
                                x = sum(x > min(x[x > 0]) * max_w_ev)
                              })

    mat_w_test = apply(X = mat_w_test, MARGIN = 2,
                       FUN = function(x){
                         x = pmin(x, min(x[x > 0]) * max_w_ev)
                         x = x / sum(x)
                       })
  }

  # Build train & test
  train = data[data$is_train == 1,]
  if (!is.null(test)){
    test = data[data$is_train == 0,]
  }


  # Calibration of RSF model
  #formula = stats::as.formula(paste0("Surv(", y_var, ",", delta_var," ) ~ ."))
  #ntime = seq(from = 0, to = max_time * 1.05, length.out = 100)

  rrtRF.fit = rrtRF(data = train,
                    y_var = y_var,
                    delta_var = delta_var,
                    x_vars = x_vars,
                    phi = phi,
                    phi.args = phi.args,
                    max_time = max_time,
                    ntree = ntree,
                    minleaf = minleaf,
                    maxdepth = maxdepth,
                    ...)

  # a corriger dans le cas de exp
  # overfitted_predictions_exp = predict_rpartRF(obj = rrtRF.fit,
  #                                          newdata = train[,x_vars],
  #                                          type = "exp")$pred

  res_overfitted_predictions_KMloc = predict_rrt_reg(obj = rrtRF.fit,
                                                     newdata = train[,x_vars],
                                                     type = "KMloc")

  overfitted_predictions_KMloc = res_overfitted_predictions_KMloc$pred_KMloc

  # Performances on train test
  # Eval on train

  # perf_train_exp = eval_model(predictions = overfitted_predictions_exp,
  #                         data = train,
  #                         phi_name = "phi",
  #                         y_name = "y_prime",
  #                         delta_name = "delta_prime",
  #                         max_time = max_time,
  #                         ev_methods = ev_methods,
  #                         phi = phi,
  #                         phi.args = phi.args,
  #                         mat_w = mat_w_train,
  #                         phi_non_censored_name = phi_non_censored_name,
  #                         bandwidths = bandwidths)


  perf_train_KMloc = eval_model(predictions = overfitted_predictions_KMloc,
                          data = train,
                          phi_name = "phi",
                          y_name = "y_prime",
                          delta_name = "delta_prime",
                          max_time = max_time,
                          ev_methods = ev_methods,
                          phi = phi,
                          phi.args = phi.args,
                          mat_w = mat_w_train,
                          phi_non_censored_name = phi_non_censored_name,
                          bandwidths = bandwidths)

  if(!is.null(test)){

    # Predictions on test set

    # test_predictions_exp = predict_rpartRF(obj = rrtRF.fit,
    #                                    newdata = test[,x_vars],
    #                                    type = "exp")$pred


    res_test_predictions_KMloc = predict_rpartRF(obj = rrtRF.fit,
                                                 newdata = test[,x_vars],
                                                 type = "KMloc")

    test_predictions_KMloc = res_test_predictions_KMloc$pred_KMloc

    # Performances on test set
    # perf_test_exp = eval_model(predictions = test_predictions_exp,
    #                              data = test,
    #                              phi_name = "phi",
    #                              y_name = "y_prime",
    #                              delta_name = "delta_prime",
    #                              max_time = max_time,
    #                              ev_methods = ev_methods,
    #                              phi = phi,
    #                              phi.args = phi.args,
    #                              mat_w = mat_w_test,
    #                              phi_non_censored_name = phi_non_censored_name,
    #                              bandwidths = bandwidths)


    perf_test_KMloc = eval_model(predictions = test_predictions_KMloc,
                           data = test,
                           phi_name = "phi",
                           y_name = "y_prime",
                           delta_name = "delta_prime",
                           max_time = max_time,
                           ev_methods = ev_methods,
                           phi = phi,
                           phi.args = phi.args,
                           mat_w = mat_w_test,
                           phi_non_censored_name = phi_non_censored_name,
                           bandwidths = bandwidths)
  }

  result = list(
    #pred_train = as.vector(overfitted_predictions_RLT),
    #perf_train = perf_train,
    perf_train_KMloc = perf_train_KMloc,
    time_points_KMloc = res_overfitted_predictions_KMloc$time,
    surv_train_KMloc = res_overfitted_predictions_KMloc$pred_surv_KMloc,
    pred_train_KMloc = overfitted_predictions_KMloc,

    train = train[,c(y_var, delta_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)],
    mat_w_train = mat_w_train,
    max_time = max_time,
    phi = phi,
    phi.args = phi.args,
    x_vars = x_vars,
    cens_rate = sum(data$delta_prime == 0) / nrow(data),
    max_w_ev = max_w_ev,
    n_w_ev_modif_train = n_w_ev_modif_train
  )
  if (!is.null(test)){
    #result$pred_test = as.vector(test_predictions_RLT)
    #result$perf_test = perf_test

    result$perf_test_KMloc = perf_test_KMloc
    result$surv_test_KMloc = res_test_predictions_KMloc$pred_surv_KMloc
    result$pred_test_KMloc = test_predictions_KMloc

    result$test = test[,c(y_var, delta_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)]
    result$mat_w_test = mat_w_test
    result$n_w_ev_modif_test = n_w_ev_modif_test
  }
  if (rrt_obj){
    result$rrt_obj = rrtRF.fit
  }
  return(result)
}



rrtRF = function(data,
                   y_var,
                   delta_var,
                   x_vars,
                   phi,
                   phi.args,
                   max_time,
                   ntree, # peut-etre a enlever plus tard
                   minleaf,
                   maxdepth,
                   ...){
  list_models = list()
  for (i in 1:ntree){
    if (i == 1){cat(paste0("tree",1,".. "))}
    if (i == 2){cat(paste0("tree",2,".. "))}
    if (i == 5){cat(paste0("tree",5,".. "))}
    if ((i %% 10) == 0){cat(paste0("tree",i,".. "))}

    sample_train = sample(x = 1:nrow(data), size = nrow(data), replace = T)
    d_train = data[sample_train, ]
    formula = stats::as.formula(paste0("Surv(", y_var, ",", delta_var," ) ~ ."))

    model = rpart::rpart(formula = formula,
                         data = d_train[, c(y_var, delta_var, x_vars)],
                         minbucket = minleaf,
                         maxdepth = maxdepth,
                         ...)

    pred_nodes = predict_nodes(obj = model, newdata = d_train[,x_vars])
    nelson_allen_estimates = do.call(rbind,
                                     args = lapply(X = unique(pred_nodes),
                                                   FUN = function(node){
                                                     fit = survival::survfit(formula = stats::as.formula(paste0("Surv(time = ", y_var,", event = ",  delta_var, ") ~ 1")),
                                                                             data = d_train[pred_nodes == node, ])
                                                     nelson_allen = cumsum(fit$n.event/fit$n.risk)
                                                     return(
                                                       c(node,
                                                         stats::approx(x = c(0,fit$time),
                                                                       y = c(0,nelson_allen),
                                                                       xout = seq(from = 0, to = max(data[,"y_prime"]) * 1.05, length.out = 100),
                                                                       method = "constant",
                                                                       rule = 2)$y)
                                                     )
                                                   }
                                     ))
    list_models[[i]] = list(model = model,
                            nelson_allen_estimates = nelson_allen_estimates,
                            time = seq(from = 0, to = max(data[,"y_prime"]) * 1.05, length.out = 100))
  }
  return(list(list_models = list_models,
              phi = phi,
              phi.args = phi.args,
              max_time = max_time
              ))
}

#' @title Compute the prediction of a model built with \code{\link{rrt_reg}}
#'
#' @description Given a model built wtih \code{\link{rrt_reg}},
#' \code{predict_rrt_reg} allows to get the predictions of the model for new
#' observations
#'
#' @param obj A list output by \code{\link{rrt_reg}}
#' @param newdata A data.frame which contains the same variables as the ones
#' used for the training
#' @param type A character string which gives the type of terminal leaf estimator to use.
#' If \cite{"exp"}, IPCW used to grow the tree are used. If \cite{"KMloc"}, a Kaplan-Meier
#' estimator is fitted in each terminal leaf
#' @return A list with the following elements :
#' \item{pred}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' (with \eqn{T' = min(T, } \code{max_time}\eqn{)}) for the observations of \code{newdata}}
#' \item{surv}{The matrix which contains the estimated values of the survival curves at
#' \code{time_points}, for the observations of \code{newdata}}
#' \item{time_points}{The vector of the time points where the survival curves
#' are evaluated}
#'
#' @seealso \code{\link{rrt_reg}}
#'
#' @export
predict_rrt_reg = function(obj, newdata, type){
  if (type == "exp"){ #nom de la methode quand on utilise le parametre theta
    preds = do.call(what = cbind, args = lapply(X = obj$list_models,
                                                FUN = function(x){rpart:::predict.rpart(obj = x$model, newdata = newdata)}))
    return(list(pred = as.vector(rowMeans(preds))))
  }
  if (type == "KMloc"){
    preds_node = do.call(what = cbind, args = lapply(X = obj$list_models,
                                                     FUN = function(x){predict_nodes(obj = x$model, newdata = newdata)}))
    ntree = ncol(preds_node)

    #not as fast as the next solution but no bug here
    #t1 = Sys.time()
    mean_nelson_allen_estimates = do.call(what = rbind,
                                          args = lapply(X = 1:dim(preds_node)[1], FUN = function(j){ # t(preds_node)
                                            mat_nelson_allen = do.call(what = rbind,
                                                                       args = lapply(X = 1:ntree, FUN = function(i){
                                                                         obj$list_models[[i]]$nelson_allen_estimates[obj$list_models[[i]]$nelson_allen_estimates[,1] == preds_node[j,i],
                                                                                                                     2:ncol(obj$list_models[[i]]$nelson_allen_estimates)]
                                                                       }))
                                            colMeans(mat_nelson_allen)
                                          }))

    pred_surv_KMloc = exp(-mean_nelson_allen_estimates)

    pred_KMloc =
      (pred_surv_KMloc[, which(obj$list_models[[1]]$time < obj$max_time)] -
         cbind(pred_surv_KMloc[, which(obj$list_models[[1]]$time < obj$max_time)[-1] ], 0)) %*%
      sapply(X = 1:length(c(obj$list_models[[1]]$time[which(obj$list_models[[1]]$time < obj$max_time)[-1]], obj$max_time)),
             FUN = function(i){do.call(obj$phi,
                                       c(list(x=c(obj$list_models[[1]]$time[which(obj$list_models[[1]]$time < obj$max_time)[-1]],
                                                  obj$max_time)[i]),
                                         obj$phi.args))})


    return(list(time = obj$list_models[[1]]$time,
                pred_surv_KMloc = pred_surv_KMloc,
                pred_KMloc = as.vector(pred_KMloc)))
  }
}



