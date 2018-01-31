
#' @title Fit classical regression model (GAM, RF) on right censored data using IPCW
#'
#' @description \code{sw_reg} is the core function of the package.
#' It implements the method we study in [Gerb. et al.] to adapt regression
#' algorithms to right censored target variable. Given a right
#' censored variable \eqn{T}, a
#' function \code{phi} and covariates \eqn{X}, \code{sw_reg}
#' aims to estimate \eqn{E[}\code{phi}\eqn{(T)|X]}. The methods is based on the
#' IPCW (Inverse probability of Censoring Weighting) principle for right
#' censored variables which is used to compensate for the censoring. Though the method
#' may generalise to many regression algorithms, \code{sw_reg}
#' only implements random forest and GAM solutions.
#' Technicaly, \code{sw_reg} is a wrapper for
#' \code{\link[randomForestSRC]{rfsrc}} and \code{\link[rpart]{rpart}}(random
#' forest) and \code{\link[mgcv]{gam}} (GAM). Different methods are available
#' to assess the quality of fit of \code{rsf_reg}.\cr \cr
#' The notations we use are :
#' \itemize{
#' \item \eqn{C} : Censoring variable
#' \item \eqn{Y = min(T, C)}
#' \item \eqn{\delta = 1_{T \le C}}  (delta)
#' }
#'
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
#' @param type_reg A character string giving the regression algorithm to use
#' in the model (default = \code{"RF"}). Other possible value is "gam".
#'
#' @param type_w A character string giving the type of IPC weights used to
#' train the regression model (default = \code{"KM"}). Other possible values are "Cox", "RSF"
#' and "unif".
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
#' @param sw_reg_obj A boolean which indicates if the random forest
#' model fitted to the training data should
#' be returned (default = \code{TRUE})
#'
#' @param cens_mod_obj A boolean which indicates if the
#' model fitted to the censoring variable to compute the weights
#' ("KM", "Cox" or "RSF") used for training should
#' be returned (default = \code{TRUE})
#'
#' @param ev_methods A vector of character strings which gives the methods that should be
#' used for the evaluation of the model (default = \code{c("concordance","weighted")}).
#' Possible choices are \code{"concordance"}, \code{"weighted"}
#' and \code{"group"}. Multiple choices are possible.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param bandwidths A vector of real numbers for the bandwidths to use for the model
#' evaluation if \code{"group"} is used
#' as an \code{ev_methods} (default = \code{NULL} : set to 50 if \code{"group"} in \code{ev_methods}).
#' Only used if \code{"group"} is used as \code{ev_methods}.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param types_w_ev A vector of character strings which gives the
#' types of weights to be used for IPCW
#' (Inverse Probability of Censoring Weighting) in the model evaluation
#' (default = \code{c("KM")} (Kaplan Meier)).
#' Possible choices are \code{"KM"}, \code{"Cox"}, \code{"RSF"} and \code{"unif"}. See
#' \emph{Details - Evaluation criteria} for more information
#'
#' @param max_w_mod A real number which gives the maximum admissible ratio
#' for the IPC weights (default = NULL : then set to \code{floor(sqrt(nrow(train))/2)})
#' used in model fitting.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param max_w_ev A real number which gives the maximum admissible ratio
#' for the IPC weights (default = 1000) used in model evaluation.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param mat_w A matrix to provide handmade IPC weights for the model
#' evaluation (default = \code{NULL}).
#' \code{mat_w} should satisfied
#' \code{nrow(mat_w) = nrow(train) + nrow(test)}
#' and a column should correspond to a type of weights
#' (multiple columns are possible).
#' Column names of \code{mat_w} may be used to specify names
#' for the provided weights (by default names will be "w1", "w2", ...)
#'
#' @param y_no_cens_var A character string which gives the name of the
#' non censored \code{y_var} (default = \code{NULL}).
#' To be used only in the context of simulated data where full about is available
#'
#' @param mode_sw_RF An integer (\code{1} or \code{2}) which specifiy the type of weighted
#' random forest to grow : \code{1} = wRF1, \code{2} = wRF2 or wRF3 (default = \code{1}).
#' See \emph{Details - Random Forest modes} for more information.
#' Only used if \code{type_reg = "RF"}
#'
#' @param ntree A positive integer which gives the number of trees to grow
#' in the forest (default = \code{100}).
#'
#' @param minleaf A positive integer indicating the minimum number of observations
#' that must be present in a (terminal) leaf (default = \code{5}).
#'
#' @param maxdepth A positive integer indicating the maximum number of layers
#' in individual trees (default = 6).
#'
#' @param mtry A positive integer indicating the number of random variables
#' to draw at each node for the split selection (default = \code{NULL}).
#' If \code{NULL}, \code{mtry} is set to \code{floor(sqrt(length(x_vars)))} by default.
#'
#' Warning (Exception to he latter statement) : if \code{mode_sw_RF = 2},
#' \code{mtry} is set to \code{length(x_vars)}
#' and can not be modified
#'
#' @param ... Additional parameter that may be pass to the regression algorithm used
#' (see \emph{Details - Additional parameters})
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
#'
#' \item \emph{Random Forest modes}
#'
#' modes correspond to different ways to take the weights into account in the random forest :
#' \itemize{
#' \item \code{mode_sw_RF = 1} : weights are computed a single time (on the
#' whole train sample) before the growing of the forest and are passed
#' to the forest as probabilities of sampling single observations for the bootstrap
#' of the random forest. This mode corresponds to \emph{wRF1} in
#' [Gerb. et al.], it internally calls the \code{\link[randomForestSRC]{rfsrc}}
#' function.
#' \item \code{mode_sw_RF = 2} : weights are computed \code{ntree} times ; for a given
#' tree, a bootsrap sample is drawn uniformly with replacement and then weights are
#' evaluated on the bootstrap sample. The tree growing procedure use then weighted error
#' as splitting criteria. Two types of predictions are made in this mode : the first
#' prediction is output as \code{pred_train/test} and it uses the same weights
#' as those used for training to compute predictions in terminal leafs (\emph{wRF2} in
#' [Gerb. et al.]). The second prediction is output as \code{pred_train/test_KMloc}
#' and it makes terminal leafs estimation by using Kaplan Meier to estimate
#' the within leaf survival function of \eqn{T}.
#' This mode internally calls the \code{\link[rpart]{rpart}} function.
#' }
#'
#' \item \emph{Additional parameters}
#'
#' \code{sw_reg} allows to pass additional parameters to
#' the underlying regression algorithm. Depending on \code{type_reg}
#' and \code{mode_sw_RF}, the wrapped function is as follow :
#' \itemize{
#' \item \code{type_reg = "RF"} and \code{mode_sw_RF = 1} : \code{\link[randomForestSRC]{rfsrc}}
#' \item \code{type_reg = "RF"} and \code{mode_sw_RF = 2} : \code{\link[rpart]{rpart}}
#' \item \code{type_reg = "gam"} : \code{\link[mgcv]{gam}}
#' }
#' For instance in the first case, one may pass to \code{sw_reg}
#' a parameter that is then passed to \code{\link[randomForestSRC]{rfsrc}}
#' (e.g. \code{proximity = TRUE})
#' }
#'
#' @return A list with the following elements :
#'
#' \item{pred_train}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' (with \eqn{T' = min(T, } \code{max_time}\eqn{)}) for the observations of the train set.
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{pred_test}{The vector of the predicted values for
#' \code{phi}\eqn{(T')} for the observations of the test set (require \code{test} != \code{NULL}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{perf_train}{The list with the values for the evaluation criteria computed on the train
#' set}
#'
#' \item{perf_test}{The list with the values for the evaluation criteria computed on the test
#' set (require \code{test} != \code{NULL}))}
#'
#' \item{w_mod_train}{The vector of the weights used to train the model,
#' after applying \code{max_w_mod} and normalising}
#'
#' \item{n_w_mod_modif_train}{The vector giving the number of train weights modified due to
#' \code{max_w_mod}}
#'
#' \item{mat_w_train}{The matrix which contains the values of the weights used for the
#' \code{"weighted"} criteria, for the observations of the train set}
#'
#' \item{mat_w_test}{The matrix which contains the values of the weights used for the
#' \code{"weighted"} criteria, for the observations of the test set}
#'
#' \item{sum_w_train}{The sum of the gross weights for the train data,
#' before applying \code{max_w_ev} and normalising}
#'
#' \item{sum_w_test}{The sum of the gross weights for the test data,
#' before applying \code{max_w_ev} and normalising}
#'
#' \item{n_w_ev_modif_train}{The vector giving the number of train weights modified due to
#' \code{max_w_ev}}
#'
#' \item{n_w_ev_modif_test}{The vector giving the number of test weights modified due to
#' \code{max_w_ev}}
#'
#' \item{sw_RF_obj}{The obj returned by \code{\link[randomForestSRC]{rfsrc}}
#' (when \code{type_reg = "RF"} and \code{mode_sw_RF = 1})}
#'
#' \item{sw_gam_obj}{The obj returned by \code{\link[mgcv]{gam}}
#' (when \code{type_reg = "gam"})}
#'
#' \item{sw_rpartRF_obj}{The obj returned by \code{\link{rfsrc}}
#' (when \code{type_reg = "RF"} and \code{mode_sw_RF = 2})}
#'
#' \item{max_time}{The real number giving the threshold used by the model}
#'
#' \item{cens_rate}{The real number giving the rate of censoring
#' of \eqn{T'}, computed on the concatenation of \code{train} and \code{test}}
#'
#' \item{pred_train_KMloc}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' for the observations of the train set (require \code{mode_sw_RF = 2}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{pred_test_KMloc}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' for the observations of the test set (require \code{mode_sw_RF = 2}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{perf_train_KMloc}{The list with the values for the evaluation criteria computed on the train
#' set (require \code{mode_sw_RF = 2}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{perf_test_KMloc}{The list with the values for the evaluation criteria computed on the test
#' set (require \code{mode_sw_RF = 2}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{surv_train_KMloc}{The matrix which contains the estimated values of the survival
#' curves at \code{time_points} (within leaf Kapaln Meier estimator),
#' for the observations of the train set (require \code{mode_sw_RF = 2})}
#'
#' \item{surv_test_KMloc}{The matrix which contains the estimated values of the survival
#' curves at \code{time_points} (within leaf Kapaln Meier estimator),
#' for the observations of the test set
#' (require \code{mode_sw_RF = 2})}
#'
#' \item{time_points}{The vector of the time points where the survival curves
#' are evaluated (require \code{mode_sw_RF = 2})}
#'
#' \item{train}{The data.frame of the train data provided as arguments, plus columns :
#' \eqn{Y' = min(Y,} \code{max_time} \eqn{)}, \eqn{\delta' = 1_{T' \le C}}
#' and \code{phi}\eqn{(T')}}
#'
#' \item{test}{The data.frame of the test data provided as arguments, plus columns :
#' \eqn{Y' = min(Y,} \code{max_time} \eqn{)}, \eqn{\delta' = 1_{T' \le C}}
#' and \code{phi}\eqn{(T')} }
#'
#' \item{type_reg}{See \emph{Argument}}
#'
#' \item{type_w}{See \emph{Argument}}
#'
#' \item{phi}{See \emph{Argument}}
#'
#' \item{phi.args}{See \emph{Argument}}
#'
#' \item{x_vars}{See \emph{Argument}}
#'
#' \item{max_w_mod}{See \emph{Argument}}
#'
#' \item{max_w_ev}{See \emph{Argument}}
#'
#' \item{mode_sw_RF}{See \emph{Argument}}
#'
#'
#' @references [Gerb. et al.] to be published
#'
#' @seealso \code{\link[randomForestSRC]{rfsrc}}, \code{\link[rpart]{rpart}},
#' \code{\link[mgcv]{gam}},
#' \code{\link{predict_sw_reg}}
#'
#' @export
#' @import methods survival
#'
#' @examples
#'
#' # ------------------------------------------------
#' #   Load "transplant" data
#' # ------------------------------------------------
#' data("transplant", package = "survival")
#' transplant$delta = 1 * (transplant$event == "ltx") # create binary var
#' # which indicate censoring/non censoring
#'
#' # keep only rows with no missing value
#' apply(transplant, MARGIN = 2, FUN = function(x){sum(is.na(x))})
#' transplant_bis = transplant[stats::complete.cases(transplant),]
#'
#' # plot the survival curve of transplant data
#' KM_transplant = survfit(formula = survival::Surv(time = futime, event = delta) ~ 1,
#'                         data = transplant_bis)
#' plot(KM_transplant)
#'
#' # ------------------------------------------------
#' #   Basic call to train a model
#' # ------------------------------------------------
#' res1 = sw_reg(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     train = transplant_bis
#' )
#' # parameters set by default
#' res1$type_w
#' res1$type_reg
#' res1$max_w_mod
#' res1$max_w_ev
#' res1$mode_sw_RF # 1 corresponds to wRF1 in [Gerb. et al.]
#'
#'
#' # train errors
#' res1$perf_train
#'
#' # ------------------------------------------------
#' #   Training with estimation of test error
#' # ------------------------------------------------
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res2 = sw_reg(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     train = transplant_bis[train_lines,],
#'                                     test = transplant_bis[-train_lines,],
#'                                     types_w_ev = c("KM", "Cox", "RSF", "unif"))
#'
#' print(res2$max_time) # default \code{max_time} has changed since train set
#' # is different
#'
#' # there is a uge overfitting in terms of quadratic errors
#' print(res2$perf_train)
#' print(res2$perf_test)
#'
#' # default parameters for the random forest are
#' res2$sw_RF_obj$ntree
#' res2$sw_RF_obj$mtry
#' res2$sw_RF_obj$nodesize
#' res2$sw_RF_obj$nodedepth # means there is no depth limit
#'
#'
#' # -----------------------------------------------------
#' #   Modify the \code{max_time} argument & look for
#' #       the best model under this setting
#' # -----------------------------------------------------
#'
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res30 = sw_reg(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     train = transplant_bis[train_lines,],
#'                                     test = transplant_bis[-train_lines,],
#'                                     type_w = "KM", # default value
#'                                     max_time = 600, # we set \code{max_time} to 600
#'                                     types_w_ev = c("KM", "Cox", "RSF", "unif"))
#' print(res30$perf_test)
#'
#' # are the other types of weights giving better results ?
#' ## Cox weights
#' res31 = sw_reg(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     train = transplant_bis[train_lines,],
#'                                     test = transplant_bis[-train_lines,],
#'                                     type_w = "Cox",
#'                                     max_time = 600, # we set \code{max_time} to 600
#'                                     types_w_ev = c("KM", "Cox", "RSF", "unif"))
#' print(res31$perf_test) # slight improvment compared with weights KM
#'
#' ## RSF weights
#' res32 = sw_reg(y_var = "futime",
#'                                      delta_var = "delta",
#'                                      x_vars = setdiff(colnames(transplant_bis),
#'                                                       c("futime", "delta", "event")),
#'                                      train = transplant_bis[train_lines,],
#'                                      test = transplant_bis[-train_lines,],
#'                                      type_w = "RSF",
#'                                      max_time = 600, # we set \code{max_time} to 600
#'                                      types_w_ev = c("KM", "Cox", "RSF", "unif"))
#' print(res32$perf_test)
#'
#' ## unif weights
#' res33 = sw_reg(y_var = "futime",
#'                                      delta_var = "delta",
#'                                      x_vars = setdiff(colnames(transplant_bis),
#'                                                       c("futime", "delta", "event")),
#'                                      train = transplant_bis[train_lines,],
#'                                      test = transplant_bis[-train_lines,],
#'                                      type_w = "unif",
#'                                      max_time = 600, # we set \code{max_time} to 600
#'                                      types_w_ev = c("KM", "Cox", "RSF", "unif"))
#' print(res33$perf_test)
#'
#' # In terms of quadratic, the best weights are the Cox weights
#' # remark : in this example there is not a big difference between the "unif" weights
#' # and the other weights because there is little censoring in the data :
#' res33$cens_rate
#'
#' # -----------------------------------------------------------------
#' #     Try wRF2  and wRF3 (both are obtained with
#' #               \code{mode_sw_RF = 2})
#' # -----------------------------------------------------------------
#'
#' res40 = sw_reg(y_var = "futime",
#'                                      delta_var = "delta",
#'                                      x_vars = setdiff(colnames(transplant_bis),
#'                                                       c("futime", "delta", "event")),
#'                                      train = transplant_bis[train_lines,],
#'                                      test = transplant_bis[-train_lines,],
#'                                      type_w = "Cox",
#'                                      max_time = 600, # we set \code{max_time} to 600
#'                                      types_w_ev = c("KM", "Cox", "RSF", "unif"),
#'                                      mode_sw_RF = 2)
#' print(res40$perf_test) # wRF2 : not as good as wRF1
#' print(res40$perf_test_KMloc) # wRF3 : worse than wRF2
#'
#'
#' # -------------------------------------------------------
#' #               Try a GAM model
#' # -------------------------------------------------------
#'
#' ## GLM with Cox weights
#' res5 = sw_reg(y_var = "futime",
#'                                      delta_var = "delta",
#'                                      x_vars = setdiff(colnames(transplant_bis),
#'                                                       c("futime", "delta", "event")),
#'                                      train = transplant_bis[train_lines,],
#'                                      test = transplant_bis[-train_lines,],
#'                                      type_w = "Cox",
#'                                      max_time = 600, # we set \code{max_time} to 600
#'                                      types_w_ev = c("KM", "Cox", "RSF", "unif"),
#'                                     type_reg = "gam")
#' print(res5$perf_test) # not as good as random forest
#'
#'
#' # ------------------------------------------------------------
#' #     Analyse the weights used for "weighted" criteria
#' # ------------------------------------------------------------
#'
#' print(res31$cens_rate) # rate of censoring taking into account \code{max_time}
#' print(head(res31$mat_w_test))
#' ## ratio max(weights)/min(weights)
#' print(apply(X = res31$mat_w_test,
#'             MARGIN = 2,
#'             FUN = function(x){max(x[x != 0])/min(x[x != 0])}))
#' # ratios are low because the censoring rate is low
#'
#' # in this case, it is not meaningful to to modify the
#' # \code{max_w_ev} argument since the maximum ratios
#' # between weights are around 2 and the test data has 197 rows.
#' # But in other situation it may be pertinent
#'
#'
#' # ------------------------------------------------
#' #   Use custom \code{phi} function
#' # ------------------------------------------------
#' g = function(x,a) abs(x-a)
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res6 = sw_reg(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     train = transplant_bis[train_lines,],
#'                                     test = transplant_bis[-train_lines,],
#'                                     phi = g,
#'                                     phi.args = list(a = 200),
#'                                     type_w = "Cox",
#'                                     max_time = 600, # we set \code{max_time} to 600
#'                                     types_w_ev = c("KM", "Cox", "RSF", "unif"))
#' print(res6$perf_test) # slight improvment compared with weights KM

sw_reg = function(y_var,
                  delta_var,
                  x_vars,
                  train,
                  test = NULL,
                  type_reg = "RF",
                  type_w = "KM",
                  phi = function(x){x},
                  phi.args = list(),
                  max_time = NULL,
                  sw_reg_obj = T,
                  cens_mod_obj = T,
                  ev_methods = c("concordance","weighted"),
                  bandwidths = NULL,
                  types_w_ev = "KM",
                  max_w_mod = NULL,
                  max_w_ev = 1000,
                  mat_w = NULL, # when mat_w is not NULL, and type_w is NULL, the column 1 of mat_w is employed for fitting
                  y_no_cens_var = NULL,

                  # param for RF
                  mode_sw_RF = 1,
                  ntree = 100,
                  minleaf = 5,
                  maxdepth = 6,
                  mtry = NULL,
                  ...){

  # Preprocessing of the arguments & data

  ## specific preprocessing compared to Cox/RSF regression
  type_w = match.arg(as.character(type_w), c("KM", "Cox", "RSF", "unif"))
  weights_manual = !is.null(mat_w)
  type_reg = match.arg(as.character(type_reg), c("RF", "gam"))
  # --------------------------------------------------------------------------------

  ev_methods <- match.arg(as.character(ev_methods), c("concordance", "group", "weighted"), several.ok = T)
  if (is.null(bandwidths) & ("group" %in% ev_methods)) bandwidths = 50
  if (is.null(max_w_mod)) max_w_mod = floor(sqrt(nrow(train))/2)

  types_w_ev = match.arg(as.character(types_w_ev), c("KM", "Cox", "RSF", "unif"), several.ok = T)
  types_w_ev = unique(c(type_w, types_w_ev)) # weights for trainig are used for evaluation

  if(is.null(mtry)){mtry = floor(sqrt(length(x_vars)))}

  # column names of mat_w should be explicit
  if (!is.null(mat_w) & is.null(colnames(mat_w))) colnames(mat_w) = paste0("w",1:ncol(mat_w))

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

  cens_mod_object = NULL

  if (weights_manual){
    if (nrow(mat_w) != nrow(data)){stop("mat_w should satisfy nrow(mat_w) = nrow(train) + nrow(test)")}
    if (is.null(type_w)){
      type_w = colnames(mat_w)[1]
    }
    w_mod_train = mat_w[1:nrow(train), type_w]
    mat_w_train = mat_w[1:nrow(train),]
    sum_w_train = apply(X = mat_w_train, MARGIN = 2, FUN = sum)
    if (!is.null(test)){
      mat_w_test = mat_w[(nrow(train) +1):nrow(data),]
      sum_w_test = apply(X = mat_w_test, MARGIN = 2, FUN = sum)
    }
  }

  if (!weights_manual){
    # Computation of the weitghts_eval
    mat_w_train = matrix(rep(0, length(types_w_ev) * sum(data$is_train == 1) ), ncol = length(types_w_ev))
    if (!is.null(test)){
      mat_w_test = matrix(rep(0, length(types_w_ev) * sum(data$is_train == 0) ), ncol = length(types_w_ev))
    }
    if (is.null(test)){
      mat_w_test = NULL
    }

    for (j in 1:length(types_w_ev)){
      if (types_w_ev[j] == type_w){
        # need to differentiate type_w or not because of the option "cens_mod_obj"
        res_weights_train = make_weights(data = data[data$is_train == 1, ],
                                         y_name = "y_prime",
                                         delta_name = "delta_prime",
                                         y_name2 = y_var,
                                         delta_name2 = delta_var,
                                         type = types_w_ev[j],
                                         max_ratio_weights = 1000,
                                         x_vars = x_vars,
                                         cens_mod_obj = cens_mod_obj)
        mat_w_train[,j] = res_weights_train$weights
        w_mod_train = res_weights_train$weights
        cens_mod_object = res_weights_train$cens_mod_object

        if (!is.null(test)){
          res_weights_test = make_weights(data = data[data$is_train == 0, ],
                                          y_name = "y_prime",
                                          delta_name = "delta_prime",
                                          y_name2 = y_var,
                                          delta_name2 = delta_var,
                                          type = types_w_ev[j],
                                          max_ratio_weights = 1000,
                                          x_vars = x_vars,
                                          cens_mod_obj = FALSE)
          mat_w_test[,j] = res_weights_test$weights
        }
      } else {
        res_weights_train = make_weights(data = data[data$is_train == 1, ],
                                         y_name = "y_prime",
                                         delta_name = "delta_prime",
                                         y_name2 = y_var,
                                         delta_name2 = delta_var,
                                         type = types_w_ev[j],
                                         max_ratio_weights = 1000,
                                         x_vars = x_vars,
                                         cens_mod_obj = FALSE)
        mat_w_train[,j] = res_weights_train$weights

        if (!is.null(test)){
          res_weights_test = make_weights(data = data[data$is_train == 0, ],
                                          y_name = "y_prime",
                                          delta_name = "delta_prime",
                                          y_name2 = y_var,
                                          delta_name2 = delta_var,
                                          type = types_w_ev[j],
                                          max_ratio_weights = 1000,
                                          x_vars = x_vars,
                                          cens_mod_obj = FALSE)
          mat_w_test[,j] = res_weights_test$weights
        }
      }
    }

    colnames(mat_w_train) = types_w_ev
    sum_w_train = apply(X = mat_w_train, MARGIN = 2, FUN = sum)
    names(sum_w_train) = types_w_ev
    if (!is.null(test)){
      colnames(mat_w_test) = types_w_ev
      sum_w_test = apply(X = mat_w_test, MARGIN = 2, FUN = sum)
      names(sum_w_test) = types_w_ev
    }
  }

  # Thresholding of the weights
  # weigths_model
  ## train
  n_w_mod_modif_train = sum(w_mod_train > (min(w_mod_train[w_mod_train > 0]) * max_w_mod))
  w_mod_train = pmin(w_mod_train, min(w_mod_train[w_mod_train > 0]) * max_w_mod)
  w_mod_train = w_mod_train / sum(w_mod_train)


  # weights_eval
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

  # build train & test
  train = data[data$is_train == 1,]
  if (!is.null(test)){
    test = data[data$is_train == 0,]
  }


  weighted_regression_result = make_res_weighted_regression(y_var = y_var,
                                                            delta_var = delta_var,
                                                            x_vars = x_vars,
                                                            train = train,
                                                            test = test,
                                                            type_reg = type_reg,
                                                            w_mod_train = w_mod_train,
                                                            phi = phi,
                                                            phi.args = phi.args,
                                                            max_time = max_time,
                                                            sw_reg_obj = sw_reg_obj,
                                                            ev_methods = ev_methods,
                                                            bandwidths = bandwidths,
                                                            mat_w_train = mat_w_train,
                                                            mat_w_test = mat_w_test,
                                                            phi_non_censored_name = phi_non_censored_name,

                                                            ## add for rpartRF
                                                            mode_sw_RF = mode_sw_RF,
                                                            type_w = type_w,
                                                            max_w_mod = max_w_mod,
                                                            ntree = ntree,
                                                            minleaf = minleaf,
                                                            maxdepth = maxdepth,
                                                            mtry = mtry,
                                                            ...)

  weighted_regression_result$type_w = type_w
  weighted_regression_result$max_w_mod = max_w_mod
  weighted_regression_result$max_w_ev = max_w_ev
  weighted_regression_result$mode_sw_RF = mode_sw_RF
  weighted_regression_result$sum_w_train = sum_w_train
  weighted_regression_result$n_w_mod_modif_train = n_w_mod_modif_train
  weighted_regression_result$n_w_ev_modif_train = n_w_ev_modif_train

  if (!is.null(test)){
    weighted_regression_result$sum_w_test = sum_w_test
    weighted_regression_result$n_w_ev_modif_test = n_w_ev_modif_test
  }
  if (!is.null(cens_mod_object)){
    # return only the model for censoring fitted on train data
    weighted_regression_result$cens_mod_obj = cens_mod_object
  }
  return(weighted_regression_result)
  }



#' @title Compute the prediction of a model built with \code{\link{sw_reg}}
#'
#' @description Given a model built wtih \code{\link{sw_reg}},
#' \code{predict_sw_reg} allows to get the
#' predictions of the model for new
#' observations
#'
#' @param obj A list output by \code{\link{sw_reg}}
#' @param newdata A data.frame which contains the same variables as the ones
#' used for the training
#' @return A list with the following elements :
#' \item{pred}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' (with \eqn{T' = min(T, } \code{max_time}\eqn{)}) for the observations of \code{newdata}}
#' \item{pred_KMloc}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' for the observations of newdata with inner Kapaln Meier weights
#' (require \code{obj} to be trained with
#' \code{mode_w_rf = 2}).
#' See \code{\link{sw_reg}} for more information}
#' \item{surv_KMloc}{The matrix which contains the estimated values of the survival curves at
#' \code{time_points} with inner Kapaln Meier weights,
#' for the observations of \code{newdata} (require \code{mode_sw_RF = 2})}
#' \item{time_points}{The vector of the time points where the survival curves
#' are evaluated (require \code{mode_sw_RF = 2})}
#'
#'
#' @seealso \code{\link{sw_reg}}
#'
#' @export
#'
#' @examples
#'
#' # ------------------------------------------------
#' #   Load "transplant" data
#' # ------------------------------------------------
#' data("transplant", package = "survival")
#' transplant$delta = 1 * (transplant$event == "ltx") # create binary var
#' # which indicate censoring/non censoring
#'
#' # keep only rows with no missing value
#' transplant_bis = transplant[stats::complete.cases(transplant),]
#'
#' # ------------------------------------------------
#' #   Basic call to train a model
#' # ------------------------------------------------
#'
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res = sw_reg(y_var = "futime",
#'                                    delta_var = "delta",
#'                                    x_vars = setdiff(colnames(transplant_bis),
#'                                                     c("futime", "delta", "event")),
#'                                    train = transplant_bis[train_lines,],
#'                                    types_w_ev = c("KM", "Cox", "RSF", "unif"),
#'                                    mode_sw_RF = 2)
#'
#' # ------------------------------------------------
#' #   Predict on new data
#' # ------------------------------------------------
#'
#' pred = predict_sw_reg(obj = res,
#'                                             newdata = transplant_bis[-train_lines,])



predict_sw_reg = function(obj, newdata){
  if (is.null(obj$sw_RF_obj) &
      is.null(obj$sw_gam_obj) &
      is.null(obj$sw_rpartRF_obj)){
    stop("to use predict_sw_reg on a sw_reg obj,
         you shoud specify sw_reg_obj = TRUE in the call of sw_reg")
  }

  if (!is.null(obj$sw_RF_obj)){
    predictions = randomForestSRC::predict.rfsrc(obj$sw_RF_obj, newdata[,obj$x_vars])$predicted
    return(list(pred = as.vector(predictions)))
  }
  if (!is.null(obj$sw_gam_obj)){
    predictions = mgcv::predict.gam(obj$sw_gam_obj, newdata[,obj$x_vars], type = "response")
    return(list(pred = as.vector(predictions)))
  }
  if (!is.null(obj$sw_rpartRF_obj)){
    res_predictions = predict_rpartRF(obj = obj$sw_rpartRF_obj,
                                      newdata = newdata[,obj$x_vars],
                                      type = "normal")

    res_predictions_KMloc = predict_rpartRF(obj = obj$sw_rpartRF_obj,
                                            newdata = newdata[,obj$x_vars],
                                            type = "KMloc")

    return(list(pred = res_predictions$pred,
                pred_KMloc = res_predictions_KMloc$pred_KMloc,
                surv_KMloc = res_predictions_KMloc$pred_surv_KMloc,
                time_points = res_predictions_KMloc$time))
  }
  }


make_res_weighted_regression = function(y_var,
                                        delta_var,
                                        x_vars,
                                        train,
                                        test,
                                        type_reg,
                                        w_mod_train,
                                        phi,
                                        phi.args,
                                        max_time,
                                        sw_reg_obj,
                                        ev_methods,
                                        bandwidths,
                                        mat_w_train,
                                        mat_w_test,
                                        phi_non_censored_name,

                                        # sup. arguments for rpartRF
                                        mode_sw_RF,
                                        type_w,
                                        max_w_mod,
                                        ntree,
                                        minleaf,
                                        maxdepth,
                                        mtry,
                                        ...){


  # Build and predict on train
  if (type_reg == "RF"){
    if (mode_sw_RF == 1){
      RF_fit = randomForestSRC::rfsrc(formula = phi ~ .,
                                      data = train[w_mod_train > 0, c("phi",x_vars)],
                                      case.wt =  w_mod_train[w_mod_train > 0],
                                      forest = T,
                                      ntree = ntree,
                                      mtry = mtry,
                                      nodesize = minleaf,
                                      nodedepth = maxdepth,
                                      ...)

      overfitted_predictions = randomForestSRC::predict.rfsrc(RF_fit,
                                                              train[,x_vars])$predicted
    }
    if (mode_sw_RF == 2){
      rpartRF_fit = rpartRF(data = train,
                            y_var = y_var,
                            delta_var = delta_var,
                            x_vars = x_vars,
                            type_w = type_w,
                            phi = phi,
                            phi.args = phi.args,
                            max_time = max_time,
                            max_w_mod = max_w_mod,
                            ntree = ntree,
                            minleaf = minleaf,
                            maxdepth = maxdepth,
                            ...)


      overfitted_predictions = predict_rpartRF(obj = rpartRF_fit,
                                               newdata = train[,x_vars],
                                               type = "normal")$pred

      res_overfitted_predictions_KMloc = predict_rpartRF(obj = rpartRF_fit,
                                                         newdata = train[,x_vars],
                                                         type = "KMloc")

      overfitted_predictions_KMloc = res_overfitted_predictions_KMloc$pred_KMloc

    }
  }

  if (type_reg == "gam"){
    gam_fit = mgcv::gam(formula = stats::as.formula(paste0("phi ~ ", paste0(x_vars, collapse = "+"))),
                        data = train[w_mod_train > 0, c("phi",x_vars)],
                        weights = w_mod_train[w_mod_train > 0],
                        ...)

    overfitted_predictions = mgcv::predict.gam(gam_fit ,
                                               train[,x_vars],
                                               type = "response")
  }

  # Eval on train
  perf_train = eval_model(predictions = overfitted_predictions,
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

  if ( (type_reg == "RF") & (mode_sw_RF == 2) ){
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
  }

  if (!is.null(test)){

    # Predict on test
    if (type_reg == "RF"){
      if (mode_sw_RF == 1){
        test_predictions = randomForestSRC::predict.rfsrc(RF_fit,
                                                          test[,x_vars])$predicted
      }
      if (mode_sw_RF == 2){

        test_predictions = predict_rpartRF(obj = rpartRF_fit,
                                           newdata = test[,x_vars],
                                           type = "normal")$pred


        res_test_predictions_KMloc = predict_rpartRF(obj = rpartRF_fit,
                                                     newdata = test[,x_vars],
                                                     type = "KMloc")

        test_predictions_KMloc = res_test_predictions_KMloc$pred_KMloc

      }
    }
    if (type_reg == "gam"){
      test_predictions = mgcv::predict.gam(gam_fit,
                                           test[,x_vars],
                                           type = "response")
    }

    # Eval on test
    perf_test = eval_model(predictions = test_predictions,
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

    if ((type_reg == "RF") & (mode_sw_RF == 2)){
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
  }

  result = list(
    pred_train = overfitted_predictions,
    perf_train = perf_train,
    train = train[,c(y_var, delta_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)],
    w_mod_train = w_mod_train,
    mat_w_train = mat_w_train,
    type_reg = type_reg,
    x_vars = x_vars,
    max_time = max_time,
    phi = phi,
    phi.args = phi.args,
    cens_rate = (sum(train$delta_prime == 0) + sum(test$delta_prime == 0)) / (nrow(train) + nrow(test))
  )
  if(sw_reg_obj){
    if (type_reg == "RF"){
      if (mode_sw_RF == 1){
        result$sw_RF_obj = RF_fit
      }
      if (mode_sw_RF == 2){
        result$sw_rpartRF_obj = rpartRF_fit
      }
    }
    if (type_reg == "gam"){
      result$sw_gam_obj = gam_fit
    }
  }
  if ((type_reg == "RF") & (mode_sw_RF == 2)){
    result$perf_train_KMloc = perf_train_KMloc
    result$time_points = res_overfitted_predictions_KMloc$time
    result$surv_train_KMloc = res_overfitted_predictions_KMloc$pred_surv_KMloc
    result$pred_train_KMloc = overfitted_predictions_KMloc
  }
  if (!is.null(test)){
    result$pred_test = test_predictions
    result$perf_test = perf_test
    result$test = test[,c(y_var, delta_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)]
    result$mat_w_test = mat_w_test
    if ((type_reg == "RF") & (mode_sw_RF == 2)){
      result$perf_test_KMloc = perf_test_KMloc
      result$surv_test_KMloc = res_test_predictions_KMloc$pred_surv_KMloc
      result$pred_test_KMloc = test_predictions_KMloc
    }
  }
  return(result)
}


predict_nodes = function (obj, newdata, na.action = stats::na.pass) {
  where <-
    if (missing(newdata))
      obj$where
  else {
    if (is.null(attr(newdata, "terms"))) {
      Terms <- stats::delete.response(obj$terms)
      newdata <- stats::model.frame(Terms, newdata, na.action = na.action,
                                    xlev = attr(obj, "xlevels"))
      if (!is.null(cl <- attr(Terms, "dataClasses")))
        stats::.checkMFClasses(cl, newdata, TRUE)
    }
    rpart:::pred.rpart(obj, rpart:::rpart.matrix(newdata))
  }
  as.integer(row.names(obj$frame))[where]
}


rpartRF = function(data,
                   y_var,
                   delta_var,
                   x_vars,
                   type_w,
                   phi,
                   phi.args,
                   max_time,
                   max_w_mod,
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

    weights_in_bag = make_weights(y_name = "y_prime",
                                  delta_name = "delta_prime",
                                  y_name2 = y_var,
                                  delta_name2 = delta_var,
                                  x_vars = x_vars,
                                  cens_mod_obj = F,
                                  max_ratio_weights = max_w_mod,
                                  type = type_w,
                                  data = d_train)$weights

    model = rpart::rpart(formula = phi ~ .,
                         data = d_train[which(weights_in_bag > 0), c("phi", x_vars)],
                         weights = weights_in_bag[which(weights_in_bag > 0)],
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
  )
  )
}

predict_rpartRF = function(obj, newdata, type){
  if (type == "normal"){
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
    # #Sys.time() - t1

    # t1 = Sys.time()
    # mean_nelson_allen_estimates = do.call(what = rbind,
    #                                       args = lapply(X = 1:dim(preds_node)[1], FUN = function(j){
    #                                         nelson_allen = rep(0, dim(obj$list_models[[1]]$nelson_allen_estimates)[2] - 1)
    #                                         for (i in 1:ntree){
    #                                           nelson_allen =
    #                                             nelson_allen +
    #                                             obj$list_models[[i]]$nelson_allen_estimates[obj$list_models[[i]]$nelson_allen_estimates[,1] == preds_node[j,i],
    #                                                                                2:ncol(obj$list_models[[i]]$nelson_allen_estimates)]
    #                                         }
    #                                         return(nelson_allen / ntree)
    #                                       }))
    # Sys.time() - t1

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
                pred_surv_KMloc = exp(-mean_nelson_allen_estimates),
                pred_KMloc = as.vector(pred_KMloc)))
  }
}


