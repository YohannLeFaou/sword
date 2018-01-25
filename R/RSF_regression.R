
#' @title Fit a RSF model and compute predictions of expectations for \code{phi}\eqn{(T)}
#'
#' @description \code{rsf_reg} is a benchmark model we use in [Gerb. et al.] (see §?).
#' To model the variable \code{phi}\eqn{(T)}, where \eqn{T} is a right censored time and
#' \code{phi} is a given function, we first
#' fit a RSF model to the data to estimate the surv function of \eqn{T} given the covariates. Then,
#' we deduce an estimator of \code{phi}\eqn{(T)} by integration of the function \code{phi} with respect
#' to the estimated surv function. Different methods are available to assess the quality of fit of
#' \code{rsf_reg}. \code{rsf_reg} is a wrapper for the \code{\link[randomForestSRC]{rfsrc}}
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
#' @param rsf_obj A boolean which indicates if the RSF model fitted to the training data should
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
#' Possible choices are \code{"KM"}, \code{"Cox"}, \code{"RSF"} and \code{"unif"}. See
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
#' @param ... Additional parameter that may be pass to the \code{\link[randomForestSRC]{rfsrc}}
#' function (package \emph{surv})
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
#' to different ways to estimate the surv function of the censoring. On the other hand, \code{"unif"}
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
#' censored variables. It indicates if the order of the pred values of the model is similar to
#' the order of the observed values. The concordance generalizes the Kendall tau to the censored case.
#'
#'
#' \item \code{"group"} : This is an experimental criteria that we didn't mention in [Gerb. et al.]. Here,
#' the idea to take the censoring into account is to measure errors given groups of observations and
#' not single observations. First, the test sample is ordered w.r.t. the pred values of the
#' model. Second, respecting this order the test sample is splited into groups of size \code{v_bandwidht[1]}, or
#' \code{v_bandwidht[2]}, etc ... (each bandwidth in \code{v_bandwidht} corresponds to a different
#' score \code{"group"}).
#'
#' Then, inside a group, an estimator of the surv function of \eqn{T} may be obtained by
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
#' \item{pred_train}{The vector of the pred values for \code{phi}\eqn{(T')}
#' (with \eqn{T' = min(T, } \code{max_time}\eqn{)}) for the observations of the train set}
#' \item{pred_test}{The vector of the pred values for
#' \code{phi}\eqn{(T')} for the observations of the test set (require \code{test} != \code{NULL})}
#' \item{perf_train}{The list with the values for the evaluation criteria computed on the train
#' set}
#' \item{perf_test}{The list with the values for the evaluation criteria computed on the test
#' set (require \code{test} != \code{NULL}))}
#' \item{surv_train}{The matrix which contains the estimated values of the surv curves at
#' \code{time_points} (with the Cox model), for the observations of the train set}
#' \item{surv_test}{The matrix which contains the estimated values of the surv curves at
#' \code{time_points}, for the observations of the test set (require \code{test} != \code{NULL})}
#' \item{time_points}{The vector of the time points where the surv curves are evaluated}
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
#' \item{rsf_obj}{The obj returned by the \code{coxph} function}
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
#' @references [Gerb. et al.] to be published
#'
#' @seealso \code{\link[randomForestSRC]{rfsrc}}, \code{\link{predict_rsf_reg}}, \url{http://rstudio.com}
#' (only here for the example)
#'
#' @export
#'
#' @examples
#'
#' # ------------------------------------------------
#' #   Load "transplant" data
#' # ------------------------------------------------
#' data("transplant", package = "surv")
#' transplant$delta = 1 * (transplant$event == "ltx") # create binary var
#' # which indicate censoring/non censoring
#'
#' # keep only rows with no missing value
#' apply(transplant, MARGIN = 2, FUN = function(x){sum(is.na(x))})
#' transplant_bis = transplant[stats::complete.cases(transplant),]
#'
#' # plot the surv curve of transplant data
#' KM_transplant = survfit(formula = surv::Surv(time = futime, event = delta) ~ 1,
#'                         data = transplant_bis)
#' plot(KM_transplant)
#'
#' # ------------------------------------------------
#' #   Basic call to train a model
#' # ------------------------------------------------
#'
#' res1 = rsf_reg(y_var = "futime",
#'                       delta_var = "delta",
#'                       x_vars = setdiff(colnames(transplant_bis),
#'                                        c("futime", "delta", "event")),
#'                       train = transplant_bis,
#'                       types_w_ev = c("KM", "Cox", "RSF", "unif"))
#'
#' # by default, (main) parameters used for the random surv forest are :
#' print(res1$rsf_obj$mtry)
#' print(res1$rsf_obj$nodesize)
#' print(res1$rsf_obj$nodedepth) # "-1" = no depth limitation
#'
#' # visualise the train predictions
#' matplot(y = t(res1$surv_train[1:30,]), x = res1$time_points, type = "l")
#' print(res1$perf_train)
#'
#' print(res1$max_time) # by default \code{max_time} is set to 2055 which is very large
#' # given the outlook of the surv function of \eqn{T}. Train errors may be
#' # overfitted
#'
#' # ------------------------------------------------
#' #   Training with estimation of test error
#' # ------------------------------------------------
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res2 = rsf_reg(y_var = "futime",
#'                       delta_var = "delta",
#'                       x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
#'                       train = transplant_bis[train_lines,],
#'                       test = transplant_bis[-train_lines,],
#'                       types_w_ev = c("KM", "Cox", "RSF", "unif"))
#'
#' print(res2$max_time) # default \code{max_time} has changed since the train set
#' # is different
#'
#' # train error is positive but test error is negative
#' print(res2$perf_train)
#' print(res2$perf_test) # weighted criterias show there is a lot of overfitting
#'
#' # ------------------------------------------------
#' #   Modify the \code{max_time} argument
#' # ------------------------------------------------
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res3 = rsf_reg(y_var = "futime",
#'                       delta_var = "delta",
#'                       x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
#'                       train = transplant_bis[train_lines,],
#'                       test = transplant_bis[-train_lines,],
#'                       max_time = 600,
#'                       types_w_ev = c("KM", "Cox", "RSF", "unif"))
#'
#' print(res3$perf_train)
#' print(res3$perf_test) # test error is much better
#'
#' # visualise the predictions
#' print(res3$pred_test[1:30])
#' matplot(y = t(res3$surv_test[1:30,]), x = res3$time_points, type = "l")
#'
#' # analyse the weights used for "weighted" criteria
#' print(res3$cens_rate) # rate of censoring taking into account \code{max_time}
#' print(head(res3$mat_w_test))
#' ## ratio max(weights)/min(weights)
#' print(apply(X = res3$mat_w_test,
#'             MARGIN = 2,
#'             FUN = function(x){max(x[x != 0])/min(x[x != 0])}))
#' # ratios are low because the censoring rate is low
#'
#' # in this case, it is not meaningful to to modify the
#' # \code{max_w_ev} argument since the maximum ratios
#' # between weights are around 2 and the test data has 197 rows.
#' # But in other situation it may be pertinent
#'
#' # ----------------------------------------------------------------
#' #   Change the parameter for \code{\link[randomForestSRC]{rfsrc}}
#' # ----------------------------------------------------------------
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res4 = rsf_reg(y_var = "futime",
#'                       delta_var = "delta",
#'                       x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
#'                       train = transplant_bis[train_lines,],
#'                       test = transplant_bis[-train_lines,],
#'                       max_time = 600,
#'                       types_w_ev = c("KM", "Cox", "RSF"),
#'                       minleaf = 5)
#'
#' # \code{minleaf} changes \code{nodesize} for the inner call to
#' # \code{\link[randomForestSRC]{rfsrc}}
#'
#' print(res4$perf_test) # slight amelioration compared
#'
#'
#' # ------------------------------------------------
#' #   Use custom \code{phi} function
#' # ------------------------------------------------
#' g = function(x,a) abs(x-a)
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res5 = rsf_reg(y_var = "futime",
#'                       delta_var = "delta",
#'                       x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
#'                       train = transplant_bis[train_lines,],
#'                       test = transplant_bis[-train_lines,],
#'                       phi = g,
#'                       phi.args = list(a = 200), # set value for "a"
#'                       max_time = 600,
#'                       types_w_ev = c("KM", "Cox", "RSF", "unif"))
#'
#' print(res5$perf_test)
#' print(res5$pred_test[1:30])

rsf_reg = function(y_var,
                          delta_var,
                          x_vars,
                          train,
                          test = NULL,
                          phi = function(x){x},
                          phi.args = list(),
                          max_time = NULL,
                          rsf_obj = T,
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

  # preprocessing() # this doesn't work for the moment

  # column names of mat_w should be explicit
  if(!is.null(mat_w) & is.null(colnames(mat_w))) colnames(mat_w) = paste0("w",1:ncol(mat_w))

  ev_methods <- match.arg(as.character(ev_methods), c("concordance", "group", "weighted"), several.ok = T)
  if (is.null(bandwidths) & ("group" %in% ev_methods)) bandwidths = 50
  types_w_ev = match.arg(as.character(types_w_ev), c("KM", "Cox", "RSF", "unif"), several.ok = T)

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
    mat_w_train = mat_w[1:nrow(train),]
    if (!is.null(test)){
      mat_w_test = mat_w[(nrow(train)+1):(nrow(data)),]
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
  formula = stats::as.formula(paste0("Surv(", y_var, ",", delta_var," ) ~ ."))
  ntime = seq(from = 0, to = max_time * 1.05, length.out = 100)

  rfSRC = randomForestSRC::rfsrc(formula = formula ,
                                 data = train[,c(y_var, delta_var, x_vars)],
                                 forest = T,
                                 ntime = ntime,
                                 ntree = ntree,
                                 nodesize = minleaf,
                                 nodedepth = maxdepth,
                                 mtry = mtry,
                                 ...)

  overfitted_predictions_direct_RSF =
    (cbind(1, rfSRC$surv[,which(rfSRC$time.interest < max_time)]) -
       cbind(rfSRC$surv[,which(rfSRC$time.interest < max_time)],0)) %*%
    sapply(X = 1:length(c(rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)),
           FUN = function(i){do.call(phi,
                                     c(list(x=c(rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)[i]),
                                       phi.args))})

  # Performances on train test
  perf_train = eval_model(predictions = overfitted_predictions_direct_RSF,
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
    test_predictions_surv_curves_direct_RSF =
      randomForestSRC::predict.rfsrc(rfSRC, test[,x_vars])$surv

    test_predictions_direct_RSF =
      ( cbind(1,test_predictions_surv_curves_direct_RSF[,which(rfSRC$time.interest < max_time)]) -
          cbind(test_predictions_surv_curves_direct_RSF[,which(rfSRC$time.interest < max_time)],0) ) %*%
      sapply(X = 1:length(c(rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)),
             FUN = function(i){do.call(phi,
                                       c(list(x=c(rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)[i]),
                                         phi.args))})

    # Performances on test set
    perf_test = eval_model(predictions = test_predictions_direct_RSF,
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
    pred_train = as.vector(overfitted_predictions_direct_RSF),
    perf_train = perf_train,
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
    result$pred_test = as.vector(test_predictions_direct_RSF)
    result$perf_test = perf_test
    result$test = test[,c(y_var, delta_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)]
    result$mat_w_test = mat_w_test
    result$n_w_ev_modif_test = n_w_ev_modif_test
  }
  if (rsf_obj){
    result$rsf_obj = rfSRC
    result$surv_train = cbind(1,rfSRC$surv[,which(rfSRC$time.interest < max_time)], 0)
    result$time_points = c(0,rfSRC$time.interest[which(rfSRC$time.interest < max_time)], max_time)
    if (!is.null(test)){
      result$surv_test =
        cbind(1,test_predictions_surv_curves_direct_RSF[,which(rfSRC$time.interest < max_time)],0)
    }
  }
  return(result)
}


#' @title Compute the prediction of a model built with \code{\link{rsf_reg}}
#'
#' @description Given a model built wtih \code{\link{rsf_reg}},
#' \code{predict_rsf_reg} allows to get the predictions of the model for new
#' observations
#'
#' @param obj A list output by \code{\link{rsf_reg}}
#' @param newdata A data.frame which contains the same variables as the ones
#' used for the training
#' @return A list with the following elements :
#' \item{pred}{The vector of the pred values for \code{phi}\eqn{(T')}
#' (with \eqn{T' = min(T, } \code{max_time}\eqn{)}) for the observations of \code{newdata}}
#' \item{surv}{The matrix which contains the estimated values of the surv curves at
#' \code{time_points}, for the observations of \code{newdata}}
#' \item{time_points}{The vector of the time points where the surv curves
#' are evaluated}
#'
#' @seealso \code{\link{rsf_reg}}
#'
#' @export
#'
#' @examples
#'
#' data("transplant", package = "surv")
#' transplant$delta = 1 * (transplant$event == "ltx") # create binary var
#' # which indicate censoring/non censoring
#'
#' # keep only rows with no missing value
#' transplant_bis = transplant[stats::complete.cases(transplant),]
#'
#'
#' # ------------------------------------------------
#' #   Basic call to train a model
#' # ------------------------------------------------
#'
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res1 = rsf_reg(y_var = "futime",
#'                       delta_var = "delta",
#'                       x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
#'                       train = transplant_bis[train_lines,],
#'                       types_w_ev = c("KM", "Cox", "RSF", "unif"))
#'
#' # ------------------------------------------------
#' #   Predict on new data
#' # ------------------------------------------------
#'
#' pred1 = predict_rsf_reg(obj = res1,
#'                        newdata = transplant_bis[-train_lines,])
#' print(pred1$pred[1:30])

predict_rsf_reg = function(obj, newdata){
  if (is.null(obj$rsf_obj)) stop("to use predict_rsf_reg on a rsf_reg obj,
                                       you shoud specify rsf_obj = TRUE in the call of rsf_reg")

  predictions_surv_curves =
    randomForestSRC::predict.rfsrc(obj$rsf_obj,
                                   newdata[,obj$x_vars])$surv

  # compute pred values for phi
  predictions =
    ( cbind(1,predictions_surv_curves[,which(obj[["rsf_obj"]][["time.interest"]] < obj$max_time)]) -
        cbind(predictions_surv_curves[,which(obj[["rsf_obj"]][["time.interest"]] < obj$max_time)],0) ) %*%
    sapply(X = 1:(sum(obj[["rsf_obj"]][["time.interest"]] < obj$max_time) + 1),
           FUN = function(i){do.call(obj$phi,
                                     c(list(x=c(obj[["rsf_obj"]][["time.interest"]][which(obj[["rsf_obj"]][["time.interest"]] < obj$max_time)],
                                                obj$max_time)[i]),
                                       obj$phi.args))})

  return(list(pred = as.vector(predictions),
              surv = cbind(1,predictions_surv_curves[,which(obj[["rsf_obj"]][["time.interest"]] < obj$max_time)],0),
              time_points = c(0,obj[["rsf_obj"]][["time.interest"]][which(obj[["rsf_obj"]][["time.interest"]] < obj$max_time)], obj$max_time)
  ))
}


