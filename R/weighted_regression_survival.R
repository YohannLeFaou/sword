
#' @title Fit classical regression model (GAM, RF) on right censored data using IPCW
#'
#' @description \code{weighted_regression_survival} is the core function of the package.
#' It implements the method we study in [Gerb. et al.] to adapt regression
#' algorithms to right censored target variable. Given a right
#' censored variable \eqn{T}, a
#' function \code{phi} and covariates \eqn{X}, \code{weighted_regression_survival}
#' aims to estimate \eqn{E[}\code{phi}\eqn{(T)|X]}. The methods is based on the
#' IPCW (Inverse probability of Censoring Weighting) principle for right
#' censored variables which is used to compensate for the censoring. Though the method
#' may generalise to many regression algorithms, \code{weighted_regression_survival}
#' only implements random forest and GAM solutions.
#' Technicaly, \code{weighted_regression_survival} is a wrapper for
#' \code{\link[randomForestSRC]{rfsrc}} and \code{\link[rpart]{rpart}}(random
#' forest) and \code{\link[mgcv]{gam}} (GAM). Different methods are available
#' to assess the quality of fit of \code{RSF_regression}.\cr \cr
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
#' @param data_train A data.frame of training observations. It must contain the column names \code{y_var},
#' \code{delta_var} and \code{x_vars}
#'
#' @param data_test A data.frame of testing observations (default = \code{NULL})
#'
#' @param type_regression A character string giving the regression algorithm to use
#' in the model (default = \code{"RF"})
#'
#' @param type_weights A character string giving the type of type of IPC weights used to
#' trian the random forest (default = \code{"KM"})
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
#' @param weighted_regression_object A boolean which indicates if the random forest
#' model fitted to the training data should
#' be returned (default = \code{TRUE})
#'
#' @param censoring_model_object A boolean which indicates if the
#' model fitted to the censoring variable to compute the weights
#' ("KM", "Cox" or "RSF") used for training should
#' be returned (default = \code{TRUE})
#'
#' @param eval_methods A vector of character strings which gives the methods that should be
#' used for the evaluation of the model (default = \code{c("concordance","weighted")}).
#' Possible choices are \code{"concordance"}, \code{"weighted"}
#' and \code{"group"}. Multiple choices are possible.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param v_bandwidth A vector of real numbers for the bandwidths to use for the model
#' evaluation if \code{"group"} is used
#' as an \code{eval_method} (default = \code{c(20)}). Only used if
#' \code{"group"} is used as \code{eval_method}.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param types_weights_eval A vector of character strings which gives the
#' types of weights to be used for IPCW
#' (Inverse Probability of Censoring Weighting) in the model evaluation
#' (default = \code{c("KM")} (Kaplan Meier)).
#' Possible choices are \code{"KM"}, \code{"Cox"}, \code{"RSF"} and \code{"unif"}. See
#' \emph{Details - Evaluation criteria} for more information
#'
#' @param max_ratio_weights_model A real number which gives the maximum admissible ratio
#' for the IPC weights (default = 1000) used in model fitting.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param max_ratio_weights_eval A real number which gives the maximum admissible ratio
#' for the IPC weights (default = 1000) used in model evaluation.
#' See \emph{Details - Evaluation criteria} for more information
#'
#' @param mat_weights A matrix to provide handmade IPC weights for the model
#' evaluation (default = \code{NULL}).
#' \code{mat_weights} should satisfied
#' \code{nrow(mat_weights) = nrow(data_train) + nrow(data_test)}
#' and a column should correspond to a type of weights
#' (multiple columns are possible).
#' Column names of \code{mat_weights} may be used to specify names
#' for the provided weights (by default names will be "w1", "w2", ...)
#'
#' @param y_non_censored_var A character string which gives the name of the
#' non censored \code{y_var} (default = \code{NULL}).
#' To be used only in the context of simulated data where full about is available
#'
#' @param mode_w_RF An integer (\code{1} or \code{2}) which specifiy the type of weighted
#' random forest to grow : \code{1} = wRF1, \code{2} = wRF2 or wRF3 (default = \code{1}).
#' See \emph{Details - Random Forest modes} for more information.
#' Only used if \code{type_regression = "RF"}
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
#' Warning (Exception to he latter statement) : if \code{mode_w_RF = 2},
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
#' The \code{types_weights_eval} argument allows the use of four kinds of IPC weights :
#' \code{"KM"}, \code{"Cox"}, \code{"RSF"} and \code{"unif"}. The first three types of weights correspond
#' to different ways to estimate the survival function of the censoring. On the other hand, \code{"unif"}
#' corresponds to \eqn{W_i = 1} for all i.
#'
#' Since the IPC weights may take too large values in some situation, \code{max_ratio_weights_eval} allows
#' to threshold the weights \eqn{(W_i)i} so that the ratio between the largest and the smallest weights doesn't
#' exceed \code{max_ratio_weights_eval}. If \eqn{W_max} is the maximum weight, considered weights are then
#' \eqn{min(W_i, W_max) / ( \sum_i min(W_i, W_max) ) }
#'
#' You can also manually provide weights to be used for IPCW with the argument \code{mat_weights} (those
#' weights will also be threshold w.r.t. \code{max_ration_weights_eval}). \code{mat_weights} should be a
#' matrix satisfying \code{nrow(mat_weights) = nrow(data_train) + nrow(data_test)}, any numbers of
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
#' \item \code{mode_w_RF = 1} : weights are computed a single time (on the
#' whole train sample) before the growing of the forest and are passed
#' to the forest as probabilities of sampling single observations for the bootstrap
#' of the random forest. This mode corresponds to \emph{wRF1} in
#' [Gerb. et al.], it internally calls the \code{\link[randomForestSRC]{rfsrc}}
#' function.
#' \item \code{mode_w_RF = 2} : weights are computed \code{ntree} times ; for a given
#' tree, a bootsrap sample is drawn uniformly with replacement and then weights are
#' evaluated on the bootstrap sample. The tree growing procedure use then weighted error
#' as splitting criteria. Two types of predictions are made in this mode : the first
#' prediction is output as \code{predicted_train/test} and it uses the same weights
#' as those used for training to compute predictions in terminal leafs (\emph{wRF2} in
#' [Gerb. et al.]). The second prediction is output as \code{predicted_train/test_KMloc}
#' and it makes terminal leafs estimation by using Kaplan Meier to estimate
#' the within leaf survival function of \eqn{T}.
#' This mode internally calls the \code{\link[rpart]{rpart}} function.
#' }
#'
#' \item \emph{Additional parameters}
#'
#' \code{weighted_regression_survival} allows to pass additional parameters to
#' the underlying regression algorithm. Depending on \code{type_regression}
#' and \code{mode_w_RF}, the wrapped function is as follow :
#' \itemize{
#' \item \code{type_regression = "RF"} and \code{mode_w_RF = 1} : \code{\link[randomForestSRC]{rfsrc}}
#' \item \code{type_regression = "RF"} and \code{mode_w_RF = 2} : \code{\link[rpart]{rpart}}
#' \item \code{type_regression = "gam"} : \code{\link[mgcv]{gam}}
#' }
#' For instance in the first case, one may pass to \code{weighted_regression_survival}
#' a parameter that is then passed to \code{\link[randomForestSRC]{rfsrc}}
#' (e.g. \code{proximity = TRUE})
#' }
#'
#' @return A list with the following elements :
#'
#' \item{predicted_train}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' (with \eqn{T' = min(T, } \code{max_time}\eqn{)}) for the observations of the train set.
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{predicted_test}{The vector of the predicted values for
#' \code{phi}\eqn{(T')} for the observations of the test set (require \code{data_test} != \code{NULL}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{list_criteria_train}{The list with the values for the evaluation criteria computed on the train
#' set}
#'
#' \item{list_criteria_test}{The list with the values for the evaluation criteria computed on the test
#' set (require \code{data_test} != \code{NULL}))}
#'
#' \item{v_weights_model_train}{The vector of the weights used to train the model,
#' after applying \code{max_ratio_weights_model} and normalising}
#'
#' \item{n_weights_model_modif_train}{The vector giving the number of train weights modified due to
#' \code{max_ratio_weights_model}}
#'
#' \item{mat_weights_train}{The matrix which contains the values of the weights used for the
#' \code{"weighted"} criteria, for the observations of the train set}
#'
#' \item{mat_weights_test}{The matrix which contains the values of the weights used for the
#' \code{"weighted"} criteria, for the observations of the test set}
#'
#' \item{sum_weights_train}{The sum of the gross weights for the train data,
#' before applying \code{max_ratio_weights_eval} and normalising}
#'
#' \item{sum_weights_test}{The sum of the gross weights for the test data,
#' before applying \code{max_ratio_weights_eval} and normalising}
#'
#' \item{n_weights_eval_modif_train}{The vector giving the number of train weights modified due to
#' \code{max_ratio_weights_eval}}
#'
#' \item{n_weights_eval_modif_test}{The vector giving the number of test weights modified due to
#' \code{max_ratio_weights_eval}}
#'
#' \item{weighted_RF_object}{The object returned by \code{\link[randomForestSRC]{rfsrc}}
#' (when \code{type_regression = "RF"} and \code{mode_w_RF = 1})}
#'
#' \item{weighted_gam_object}{The object returned by \code{\link[mgcv]{gam}}
#' (when \code{type_regression = "gam"})}
#'
#' \item{weighted_rpartRF_object}{The object returned by \code{\link{rfsrc}}
#' (when \code{type_regression = "RF"} and \code{mode_w_RF = 2})}
#'
#' \item{max_time}{The real number giving the threshold used by the model}
#'
#' \item{censoring_rate_with_threshold}{The real number giving the rate of censoring
#' of \eqn{T'}, computed on the concatenation of \code{data_train} and \code{data_test}}
#'
#' \item{predicted_train_KMloc}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' for the observations of the train set (require \code{mode_w_RF = 2}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{predicted_test_KMloc}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' for the observations of the test set (require \code{mode_w_RF = 2}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{list_criteria_train_KMloc}{The list with the values for the evaluation criteria computed on the train
#' set (require \code{mode_w_RF = 2}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{list_criteria_test_KMloc}{The list with the values for the evaluation criteria computed on the test
#' set (require \code{mode_w_RF = 2}).
#' See \emph{Details - Random Forest modes} for more information}
#'
#' \item{survival_train_KMloc}{The matrix which contains the estimated values of the survival
#' curves at \code{time_points} (within leaf Kapaln Meier estimator),
#' for the observations of the train set (require \code{mode_w_RF = 2})}
#'
#' \item{survival_test_KMloc}{The matrix which contains the estimated values of the survival
#' curves at \code{time_points} (within leaf Kapaln Meier estimator),
#' for the observations of the test set
#' (require \code{mode_w_RF = 2})}
#'
#' \item{time_points}{The vector of the time points where the survival curves
#' are evaluated (require \code{mode_w_RF = 2})}
#'
#' \item{data_train}{The data.frame of the train data provided as arguments, plus columns :
#' \eqn{Y' = min(Y,} \code{max_time} \eqn{)}, \eqn{\delta' = 1_{T' \le C}}
#' and \code{phi}\eqn{(T')}}
#'
#' \item{data_test}{The data.frame of the test data provided as arguments, plus columns :
#' \eqn{Y' = min(Y,} \code{max_time} \eqn{)}, \eqn{\delta' = 1_{T' \le C}}
#' and \code{phi}\eqn{(T')} }
#'
#' \item{type_regression}{See \emph{Argument}}
#'
#' \item{type_weights}{See \emph{Argument}}
#'
#' \item{phi}{See \emph{Argument}}
#'
#' \item{phi.args}{See \emph{Argument}}
#'
#' \item{x_vars}{See \emph{Argument}}
#'
#' \item{max_ratio_weights_model}{See \emph{Argument}}
#'
#' \item{max_ratio_weights_eval}{See \emph{Argument}}
#'
#' \item{mode_w_RF}{See \emph{Argument}}
#'
#'
#' @references [Gerb. et al.] to be published
#'
#' @seealso \code{\link[randomForestSRC]{rfsrc}}, \code{\link[rpart]{rpart}},
#' \code{\link[mgcv]{gam}},
#' \code{\link{predict_weighted_regression_survival}}
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
#' res1 = weighted_regression_survival(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     data_train = transplant_bis
#' )
#' # parameters set by default
#' res1$type_weights
#' res1$type_regression
#' res1$max_ratio_weights_model
#' res1$max_ratio_weights_eval
#' res1$mode_w_RF # 1 corresponds to wRF1 in [Gerb. et al.]
#'
#'
#' # train errors
#' res1$list_criteria_train
#'
#' # ------------------------------------------------
#' #   Training with estimation of test error
#' # ------------------------------------------------
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res2 = weighted_regression_survival(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     data_train = transplant_bis[train_lines,],
#'                                     data_test = transplant_bis[-train_lines,],
#'                                     types_weights_eval = c("KM", "Cox", "RSF", "unif"))
#'
#' print(res2$max_time) # default \code{max_time} has changed since train set
#' # is different
#'
#' # there is a uge overfitting in terms of quadratic errors
#' print(res2$list_criteria_train)
#' print(res2$list_criteria_test)
#'
#' # default parameters for the random forest are
#' res2$weighted_RF_object$ntree
#' res2$weighted_RF_object$mtry
#' res2$weighted_RF_object$nodesize
#' res2$weighted_RF_object$nodedepth # means there is no depth limit
#'
#'
#' # -----------------------------------------------------
#' #   Modify the \code{max_time} argument & look for
#' #       the best model under this setting
#' # -----------------------------------------------------
#'
#' set.seed(17)
#' train_lines = sample(1:nrow(transplant_bis), 600)
#' res30 = weighted_regression_survival(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     data_train = transplant_bis[train_lines,],
#'                                     data_test = transplant_bis[-train_lines,],
#'                                     type_weights = "KM", # default value
#'                                     max_time = 600, # we set \code{max_time} to 600
#'                                     types_weights_eval = c("KM", "Cox", "RSF", "unif"))
#' print(res30$list_criteria_test)
#'
#' # are the other types of weights giving better results ?
#' ## Cox weights
#' res31 = weighted_regression_survival(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     data_train = transplant_bis[train_lines,],
#'                                     data_test = transplant_bis[-train_lines,],
#'                                     type_weights = "Cox",
#'                                     max_time = 600, # we set \code{max_time} to 600
#'                                     types_weights_eval = c("KM", "Cox", "RSF", "unif"))
#' print(res31$list_criteria_test) # slight improvment compared with weights KM
#'
#' ## RSF weights
#' res32 = weighted_regression_survival(y_var = "futime",
#'                                      delta_var = "delta",
#'                                      x_vars = setdiff(colnames(transplant_bis),
#'                                                       c("futime", "delta", "event")),
#'                                      data_train = transplant_bis[train_lines,],
#'                                      data_test = transplant_bis[-train_lines,],
#'                                      type_weights = "RSF",
#'                                      max_time = 600, # we set \code{max_time} to 600
#'                                      types_weights_eval = c("KM", "Cox", "RSF", "unif"))
#' print(res32$list_criteria_test)
#'
#' ## unif weights
#' res33 = weighted_regression_survival(y_var = "futime",
#'                                      delta_var = "delta",
#'                                      x_vars = setdiff(colnames(transplant_bis),
#'                                                       c("futime", "delta", "event")),
#'                                      data_train = transplant_bis[train_lines,],
#'                                      data_test = transplant_bis[-train_lines,],
#'                                      type_weights = "unif",
#'                                      max_time = 600, # we set \code{max_time} to 600
#'                                      types_weights_eval = c("KM", "Cox", "RSF", "unif"))
#' print(res33$list_criteria_test)
#'
#' # In terms of quadratic, the best weights are the Cox weights
#' # remark : in this example there is not a big difference between the "unif" weights
#' # and the other weights because there is little censoring in the data :
#' res33$censoring_rate_with_threshold
#'
#' # -----------------------------------------------------------------
#' #     Try wRF2  and wRF3 (both are obtained with
#' #               \code{mode_w_RF = 2})
#' # -----------------------------------------------------------------
#'
#' res40 = weighted_regression_survival(y_var = "futime",
#'                                      delta_var = "delta",
#'                                      x_vars = setdiff(colnames(transplant_bis),
#'                                                       c("futime", "delta", "event")),
#'                                      data_train = transplant_bis[train_lines,],
#'                                      data_test = transplant_bis[-train_lines,],
#'                                      type_weights = "Cox",
#'                                      max_time = 600, # we set \code{max_time} to 600
#'                                      types_weights_eval = c("KM", "Cox", "RSF", "unif"),
#'                                      mode_w_RF = 2)
#' print(res40$list_criteria_test) # wRF2 : not as good as wRF1
#' print(res40$list_criteria_test_KMloc) # wRF3 : worse than wRF2
#'
#'
#' # -------------------------------------------------------
#' #               Try a GAM model
#' # -------------------------------------------------------
#'
#' ## GLM with Cox weights
#' res5 = weighted_regression_survival(y_var = "futime",
#'                                      delta_var = "delta",
#'                                      x_vars = setdiff(colnames(transplant_bis),
#'                                                       c("futime", "delta", "event")),
#'                                      data_train = transplant_bis[train_lines,],
#'                                      data_test = transplant_bis[-train_lines,],
#'                                      type_weights = "Cox",
#'                                      max_time = 600, # we set \code{max_time} to 600
#'                                      types_weights_eval = c("KM", "Cox", "RSF", "unif"),
#'                                     type_regression = "gam")
#' print(res5$list_criteria_test) # not as good as random forest
#'
#'
#' # ------------------------------------------------------------
#' #     Analyse the weights used for "weighted" criteria
#' # ------------------------------------------------------------
#'
#' print(res31$censoring_rate_with_threshold) # rate of censoring taking into account \code{max_time}
#' print(head(res31$mat_weights_test))
#' ## ratio max(weights)/min(weights)
#' print(apply(X = res31$mat_weights_test,
#'             MARGIN = 2,
#'             FUN = function(x){max(x[x != 0])/min(x[x != 0])}))
#' # ratios are low because the censoring rate is low
#'
#' # in this case, it is not meaningful to to modify the
#' # \code{max_ratio_weights_eval} argument since the maximum ratios
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
#' res6 = weighted_regression_survival(y_var = "futime",
#'                                     delta_var = "delta",
#'                                     x_vars = setdiff(colnames(transplant_bis),
#'                                                      c("futime", "delta", "event")),
#'                                     data_train = transplant_bis[train_lines,],
#'                                     data_test = transplant_bis[-train_lines,],
#'                                     phi = g,
#'                                     phi.args = list(a = 200),
#'                                     type_weights = "Cox",
#'                                     max_time = 600, # we set \code{max_time} to 600
#'                                     types_weights_eval = c("KM", "Cox", "RSF", "unif"))
#' print(res6$list_criteria_test) # slight improvment compared with weights KM

weighted_regression_survival = function(y_var,
                                        delta_var,
                                        x_vars,
                                        data_train,
                                        data_test = NULL,
                                        type_regression = "RF",
                                        type_weights = "KM",
                                        phi = function(x){x},
                                        phi.args = list(),
                                        max_time = NULL,
                                        weighted_regression_object = T,
                                        censoring_model_object = T,
                                        eval_methods = c("concordance","weighted"),
                                        v_bandwidth = 20,
                                        types_weights_eval = "KM",
                                        max_ratio_weights_model = 20,
                                        max_ratio_weights_eval = 1000,
                                        mat_weights = NULL, # we say that when mat_weights is not NULL, and type_weights is NULL, the column 1 of mat_weights is employed for fitting
                                        y_non_censored_var = NULL,

                                        # param for RF
                                        mode_w_RF = 1,
                                        ntree = 100,
                                        minleaf = 5,
                                        maxdepth = 6,
                                        mtry = NULL,
                                        ...){

  # Preprocessing of the arguments & data

  ## specific preprocessing compared to Cox/RSF regression
  type_weights = match.arg(as.character(type_weights), c("KM", "Cox", "RSF", "unif"))
  weights_manual = !is.null(mat_weights)
  type_regression = match.arg(as.character(type_regression), c("RF", "gam"))
  # --------------------------------------------------------------------------------

  eval_methods <- match.arg(as.character(eval_methods), c("concordance","single", "group", "weighted"), several.ok = T)
  types_weights_eval = match.arg(as.character(types_weights_eval), c("KM", "Cox", "RSF", "unif"), several.ok = T)
  types_weights_eval = unique(c(type_weights, types_weights_eval)) # weights for trainig are used for evaluation

  if(is.null(mtry)){mtry = floor(sqrt(length(x_vars)))}

  # column names of mat_weights should be explicit
  if (!is.null(mat_weights) & is.null(colnames(mat_weights))) colnames(mat_weights) = paste0("w",1:ncol(mat_weights))

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

  if (weights_manual){
    if (nrow(mat_weights) != nrow(data)){stop("mat_weights should satisfy nrow(mat_weights) = nrow(data_train) + nrow(data_test)")}
    if (is.null(type_weights)){
      type_weights = colnames(mat_weights)[1]
    }
    v_weights_model_train = mat_weights[1:nrow(data_train), type_weights]
    mat_weights_train = mat_weights[1:nrow(data_train),]
    sum_weights_train = apply(X = mat_weights_train, MARGIN = 2, FUN = sum)
    censoring_model_object = NULL
    if (!is.null(data_test)){
      mat_weights_test = mat_weights[(nrow(data_train) +1):nrow(data),]
      sum_weights_test = apply(X = mat_weights_test, MARGIN = 2, FUN = sum)
    }
  }

  if (!weights_manual){
    # Computation of the weitghts_eval
    mat_weights_train = matrix(rep(0, length(types_weights_eval) * sum(data$is_train == 1) ), ncol = length(types_weights_eval))
    if (!is.null(data_test)){
      mat_weights_test = matrix(rep(0, length(types_weights_eval) * sum(data$is_train == 0) ), ncol = length(types_weights_eval))
    }
    if (is.null(data_test)){
      mat_weights_test = NULL
    }

    for (j in 1:length(types_weights_eval)){
      if (types_weights_eval[j] == type_weights){
        # need to differentiate type_weights or not because of the option "censoring_model_object"
        res_weights_train = make_weights(data = data[data$is_train == 1, ],
                                         y_name = "y_prime",
                                         delta_name = "delta_prime",
                                         y_name2 = y_var,
                                         delta_name2 = delta_var,
                                         type = types_weights_eval[j],
                                         max_ratio_weights = 1000,
                                         x_vars = x_vars,
                                         censoring_model_object = censoring_model_object)
        mat_weights_train[,j] = res_weights_train$weights
        v_weights_model_train = res_weights_train$weights
        censoring_model_object = res_weights_train$censoring_model_object

        if (!is.null(data_test)){
          res_weights_test = make_weights(data = data[data$is_train == 0, ],
                                          y_name = "y_prime",
                                          delta_name = "delta_prime",
                                          y_name2 = y_var,
                                          delta_name2 = delta_var,
                                          type = types_weights_eval[j],
                                          max_ratio_weights = 1000,
                                          x_vars = x_vars,
                                          censoring_model_object = FALSE)
          mat_weights_test[,j] = res_weights_test$weights
        }
      } else {
        res_weights_train = make_weights(data = data[data$is_train == 1, ],
                                         y_name = "y_prime",
                                         delta_name = "delta_prime",
                                         y_name2 = y_var,
                                         delta_name2 = delta_var,
                                         type = types_weights_eval[j],
                                         max_ratio_weights = 1000,
                                         x_vars = x_vars,
                                         censoring_model_object = FALSE)
        mat_weights_train[,j] = res_weights_train$weights

        if (!is.null(data_test)){
          res_weights_test = make_weights(data = data[data$is_train == 0, ],
                                          y_name = "y_prime",
                                          delta_name = "delta_prime",
                                          y_name2 = y_var,
                                          delta_name2 = delta_var,
                                          type = types_weights_eval[j],
                                          max_ratio_weights = 1000,
                                          x_vars = x_vars,
                                          censoring_model_object = FALSE)
          mat_weights_test[,j] = res_weights_test$weights
        }
      }
    }

    colnames(mat_weights_train) = types_weights_eval
    sum_weights_train = apply(X = mat_weights_train, MARGIN = 2, FUN = sum)
    names(sum_weights_train) = types_weights_eval
    if (!is.null(data_test)){
      colnames(mat_weights_test) = types_weights_eval
      sum_weights_test = apply(X = mat_weights_test, MARGIN = 2, FUN = sum)
      names(sum_weights_test) = types_weights_eval
    }
  }

  # Thresholding of the weights
  # weigths_model
  ## train
  n_weights_model_modif_train = sum(v_weights_model_train > (min(v_weights_model_train[v_weights_model_train > 0]) * max_ratio_weights_model))
  v_weights_model_train = pmin(v_weights_model_train, min(v_weights_model_train[v_weights_model_train > 0]) * max_ratio_weights_model)
  v_weights_model_train = v_weights_model_train / sum(v_weights_model_train)


  # weights_eval
  ## train
  n_weights_eval_modif_train = apply(X = mat_weights_train, MARGIN = 2,
                                     FUN = function(x){
                                       x = sum(x > min(x[x > 0]) * max_ratio_weights_eval)
                                     })

  mat_weights_train = apply(X = mat_weights_train, MARGIN = 2,
                            FUN = function(x){
                              x = pmin(x, min(x[x > 0]) * max_ratio_weights_eval)
                              x = x / sum(x)
                            })
  ## test
  if (!is.null(data_test)){

    n_weights_eval_modif_test = apply(X = mat_weights_test, MARGIN = 2,
                                       FUN = function(x){
                                         x = sum(x > min(x[x > 0]) * max_ratio_weights_eval)
                                       })

    mat_weights_test = apply(X = mat_weights_test, MARGIN = 2,
                             FUN = function(x){
                               x = pmin(x, min(x[x > 0]) * max_ratio_weights_eval)
                               x = x / sum(x)
                             })
  }

  # build train & test
  data_train = data[data$is_train == 1,]
  if (!is.null(data_test)){
    data_test = data[data$is_train == 0,]
  }


  weighted_regression_result = make_res_weighted_regression(y_var = y_var,
                                                            delta_var = delta_var,
                                                            x_vars = x_vars,
                                                            data_train = data_train,
                                                            data_test = data_test,
                                                            type_regression = type_regression,
                                                            v_weights_model_train = v_weights_model_train,
                                                            phi = phi,
                                                            phi.args = phi.args,
                                                            max_time = max_time,
                                                            weighted_regression_object = weighted_regression_object,
                                                            eval_methods = eval_methods,
                                                            v_bandwidth = v_bandwidth,
                                                            mat_weights_train = mat_weights_train,
                                                            mat_weights_test = mat_weights_test,
                                                            phi_non_censored_name = phi_non_censored_name,

                                                            ## add for rpartRF
                                                            mode_w_RF = mode_w_RF,
                                                            type_weights = type_weights,
                                                            max_ratio_weights_model = max_ratio_weights_model,
                                                            ntree = ntree,
                                                            minleaf = minleaf,
                                                            maxdepth = maxdepth,
                                                            mtry = mtry,
                                                            ...)

  weighted_regression_result$type_weights = type_weights
  weighted_regression_result$max_ratio_weights_model = max_ratio_weights_model
  weighted_regression_result$max_ratio_weights_eval = max_ratio_weights_eval
  weighted_regression_result$mode_w_RF = mode_w_RF
  weighted_regression_result$sum_weights_train = sum_weights_train
  weighted_regression_result$n_weights_model_modif_train = n_weights_model_modif_train
  weighted_regression_result$n_weights_eval_modif_train = n_weights_eval_modif_train

  if (!is.null(data_test)){
    weighted_regression_result$sum_weights_test = sum_weights_test
    weighted_regression_result$n_weights_eval_modif_test = n_weights_eval_modif_test
  }
  if (!is.null(censoring_model_object)){
    # return only the model for censoring fitted on train data
    weighted_regression_result$censoring_model_object = censoring_model_object
  }
  return(weighted_regression_result)
}



#' @title Compute the prediction of a model built with \code{\link{weighted_regression_survival}}
#'
#' @description Given a model built wtih \code{\link{weighted_regression_survival}},
#' \code{predict_weighted_regression_survival} allows to get the
#' predictions of the model for new
#' observations
#'
#' @param object A list output by \code{\link{weighted_regression_survival}}
#' @param newdata A data.frame which contains the same variables as the ones
#' used for the training
#' @return A list with the following elements :
#' \item{predicted}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' (with \eqn{T' = min(T, } \code{max_time}\eqn{)}) for the observations of \code{newdata}}
#' \item{predicted_KMloc}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' for the observations of newdata with inner Kapaln Meier weights
#' (require \code{object} to be trained with
#' \code{mode_w_rf = 2}).
#' See \code{\link{weighted_regression_survival}} for more information}
#' \item{survival_KMloc}{The matrix which contains the estimated values of the survival curves at
#' \code{time_points} with inner Kapaln Meier weights,
#' for the observations of \code{newdata} (require \code{mode_w_RF = 2})}
#' \item{time_points}{The vector of the time points where the survival curves
#' are evaluated (require \code{mode_w_RF = 2})}
#'
#'
#' @seealso \code{\link{weighted_regression_survival}}
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
#' res = weighted_regression_survival(y_var = "futime",
#'                                    delta_var = "delta",
#'                                    x_vars = setdiff(colnames(transplant_bis),
#'                                                     c("futime", "delta", "event")),
#'                                    data_train = transplant_bis[train_lines,],
#'                                    types_weights_eval = c("KM", "Cox", "RSF", "unif"),
#'                                    mode_w_RF = 2)
#'
#' # ------------------------------------------------
#' #   Predict on new data
#' # ------------------------------------------------
#'
#' pred = predict_weighted_regression_survival(object = res,
#'                                             newdata = transplant_bis[-train_lines,])



predict_weighted_regression_survival = function(object, newdata){
  if (is.null(object$weighted_RF_object) &
      is.null(object$weighted_gam_object) &
      is.null(object$weighted_rpartRF_object)){
    stop("to use predict_weighted_regression_survival on a weighted_regression_survival object,
         you shoud specify weighted_regression_object = TRUE in the call of weighted_regression_survival")
  }

  if (!is.null(object$weighted_RF_object)){
    predictions = randomForestSRC::predict.rfsrc(object$weighted_RF_object, newdata[,object$x_vars])$predicted
    return(list(predicted = as.vector(predictions)))
  }
  if (!is.null(object$weighted_gam_object)){
    predictions = mgcv::predict.gam(object$weighted_gam_object, newdata[,object$x_vars])
    return(list(predicted = as.vector(predictions)))
  }
  if (!is.null(object$weighted_rpartRF_object)){
    res_predictions = predict_rpartRF(object = object$weighted_rpartRF_object,
                                  newdata = newdata[,object$x_vars],
                                  type = "normal")

    res_predictions_KMloc = predict_rpartRF(object = object$weighted_rpartRF_object,
                                       newdata = newdata[,object$x_vars],
                                       type = "KMloc")

    return(list(predicted = res_predictions$predicted,
                predicted_KMloc = res_predictions_KMloc$predicted_KMloc,
                survival_KMloc = res_predictions_KMloc$preds_surv_KMloc,
                time_points = res_predictions_KMloc$time))
  }
}


make_res_weighted_regression = function(y_var,
                                        delta_var,
                                        x_vars,
                                        data_train,
                                        data_test,
                                        type_regression,
                                        v_weights_model_train,
                                        phi,
                                        phi.args,
                                        max_time,
                                        weighted_regression_object,
                                        eval_methods,
                                        v_bandwidth,
                                        mat_weights_train,
                                        mat_weights_test,
                                        phi_non_censored_name,

                                        # sup. arguments for rpartRF
                                        mode_w_RF,
                                        type_weights,
                                        max_ratio_weights_model,
                                        ntree,
                                        minleaf,
                                        maxdepth,
                                        mtry,
                                        ...){


  # Build and predict on train
  if (type_regression == "RF"){
    if (mode_w_RF == 1){
      RF_fit = randomForestSRC::rfsrc(formula = phi ~ .,
                                      data = data_train[v_weights_model_train > 0, c("phi",x_vars)],
                                      case.wt =  v_weights_model_train[v_weights_model_train > 0],
                                      forest = T,
                                      ntree = ntree,
                                      mtry = mtry,
                                      nodesize = minleaf,
                                      nodedepth = maxdepth,
                                      ...)

      overfitted_predictions = randomForestSRC::predict.rfsrc(RF_fit,
                                                              data_train[,x_vars])$predicted
    }
    if (mode_w_RF == 2){
      rpartRF_fit = rpartRF(data = data_train,
                            y_var = y_var,
                            delta_var = delta_var,
                            x_vars = x_vars,
                            type_weights = type_weights,
                            phi = phi,
                            phi.args = phi.args,
                            max_time = max_time,
                            max_ratio_weights_model = max_ratio_weights_model,
                            ntree = ntree,
                            minleaf = minleaf,
                            maxdepth = maxdepth,
                            ...)


      overfitted_predictions = predict_rpartRF(object = rpartRF_fit,
                                               newdata = data_train[,x_vars],
                                               type = "normal")$predicted

      res_overfitted_predictions_KMloc = predict_rpartRF(object = rpartRF_fit,
                                                       newdata = data_train[,x_vars],
                                                       type = "KMloc")

      overfitted_predictions_KMloc = res_overfitted_predictions_KMloc$predicted_KMloc

    }
  }

  if (type_regression == "gam"){
    gam_fit = mgcv::gam(formula = stats::as.formula(paste0("phi ~ ", paste0(x_vars, collapse = "+"))),
                        data = data_train[v_weights_model_train > 0, c("phi",x_vars)],
                        weights = v_weights_model_train[v_weights_model_train > 0],
                        ...)

    overfitted_predictions = mgcv::predict.gam(gam_fit ,
                                               data_train[,x_vars])
  }

  # Eval on train
  list_criteria_train = eval_model(predictions = overfitted_predictions,
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

  if ( (type_regression == "RF") & (mode_w_RF == 2) ){
    list_criteria_train_KMloc = eval_model(predictions = overfitted_predictions_KMloc,
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
  }

  if (!is.null(data_test)){

    # Predict on test
    if (type_regression == "RF"){
      if (mode_w_RF == 1){
        test_predictions = randomForestSRC::predict.rfsrc(RF_fit,
                                                          data_test[,x_vars])$predicted
      }
      if (mode_w_RF == 2){

        test_predictions = predict_rpartRF(object = rpartRF_fit,
                                           newdata = data_test[,x_vars],
                                           type = "normal")$predicted


        res_test_predictions_KMloc = predict_rpartRF(object = rpartRF_fit,
                                                   newdata = data_test[,x_vars],
                                                   type = "KMloc")

        test_predictions_KMloc = res_test_predictions_KMloc$predicted_KMloc

      }
    }
    if (type_regression == "gam"){
      test_predictions = mgcv::predict.gam(gam_fit,
                                           data_test[,x_vars])
    }

    # Eval on test
    list_criteria_test = eval_model(predictions = test_predictions,
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

    if ((type_regression == "RF") & (mode_w_RF == 2)){
      list_criteria_test_KMloc = eval_model(predictions = test_predictions_KMloc,
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
  }

  result = list(
    predicted_train = overfitted_predictions,
    list_criteria_train = list_criteria_train,
    data_train = data_train[,c(y_var, delta_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)],
    v_weights_model_train = v_weights_model_train,
    mat_weights_train = mat_weights_train,
    type_regression = type_regression,
    x_vars = x_vars,
    max_time = max_time,
    phi = phi,
    phi.args = phi.args,
    censoring_rate_with_threshold = (sum(data_train$delta_prime == 0) + sum(data_test$delta_prime == 0)) / (nrow(data_train) + nrow(data_test))
  )
  if(weighted_regression_object){
    if (type_regression == "RF"){
      if (mode_w_RF == 1){
        result$weighted_RF_object = RF_fit
      }
      if (mode_w_RF == 2){
        result$weighted_rpartRF_object = rpartRF_fit
      }
    }
    if (type_regression == "gam"){
      result$weighted_gam_object = gam_fit
    }
  }
  if ((type_regression == "RF") & (mode_w_RF == 2)){
    result$list_criteria_train_KMloc = list_criteria_train_KMloc
    result$time_points = res_overfitted_predictions_KMloc$time
    result$survival_train_KMloc = res_overfitted_predictions_KMloc$preds_surv_KMloc
    result$predicted_train_KMloc = overfitted_predictions_KMloc
  }
  if (!is.null(data_test)){
    result$predicted_test = test_predictions
    result$list_criteria_test = list_criteria_test
    result$data_test = data_test[,c(y_var, delta_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)]
    result$mat_weights_test = mat_weights_test
    if ((type_regression == "RF") & (mode_w_RF == 2)){
      result$list_criteria_test_KMloc = list_criteria_test_KMloc
      result$survival_test_KMloc = res_test_predictions_KMloc$preds_surv_KMloc
      result$predicted_test_KMloc = test_predictions_KMloc
    }
  }
  return(result)
}


predict_nodes = function (object, newdata, na.action = stats::na.pass) {
  where <-
    if (missing(newdata))
      object$where
  else {
    if (is.null(attr(newdata, "terms"))) {
      Terms <- stats::delete.response(object$terms)
      newdata <- stats::model.frame(Terms, newdata, na.action = na.action,
                             xlev = attr(object, "xlevels"))
      if (!is.null(cl <- attr(Terms, "dataClasses")))
        stats::.checkMFClasses(cl, newdata, TRUE)
    }
    rpart:::pred.rpart(object, rpart:::rpart.matrix(newdata))
  }
  as.integer(row.names(object$frame))[where]
}


rpartRF = function(data,
                   y_var,
                   delta_var,
                   x_vars,
                   type_weights,
                   phi,
                   phi.args,
                   max_time,
                   max_ratio_weights_model,
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
                                  censoring_model_object = F,
                                  max_ratio_weights = max_ratio_weights_model,
                                  type = type_weights,
                                  data = d_train)$weights

    model = rpart::rpart(formula = phi ~ .,
                         data = d_train[which(weights_in_bag > 0), c("phi", x_vars)],
                         weights = weights_in_bag[which(weights_in_bag > 0)],
                         minbucket = minleaf,
                         maxdepth = maxdepth,
                         ...)

    pred_nodes = predict_nodes(object = model, newdata = d_train[,x_vars])
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

predict_rpartRF = function(object, newdata, type){
  if (type == "normal"){
    preds = do.call(what = cbind, args = lapply(X = object$list_models,
                                                FUN = function(x){rpart:::predict.rpart(object = x$model, newdata = newdata)}))
    return(list(predicted = as.vector(rowMeans(preds))))
  }
  if (type == "KMloc"){
    preds_node = do.call(what = cbind, args = lapply(X = object$list_models,
                                                     FUN = function(x){predict_nodes(object = x$model, newdata = newdata)}))
    ntree = ncol(preds_node)

    #not as fast as the next solution but no bug here
    #t1 = Sys.time()
    mean_nelson_allen_estimates = do.call(what = rbind,
            args = lapply(X = 1:dim(preds_node)[1], FUN = function(j){ # t(preds_node)
              mat_nelson_allen = do.call(what = rbind,
                      args = lapply(X = 1:ntree, FUN = function(i){
                        object$list_models[[i]]$nelson_allen_estimates[object$list_models[[i]]$nelson_allen_estimates[,1] == preds_node[j,i],
                                                           2:ncol(object$list_models[[i]]$nelson_allen_estimates)]
                      }))
              colMeans(mat_nelson_allen)
    }))
    # #Sys.time() - t1

    # t1 = Sys.time()
    # mean_nelson_allen_estimates = do.call(what = rbind,
    #                                       args = lapply(X = 1:dim(preds_node)[1], FUN = function(j){
    #                                         nelson_allen = rep(0, dim(object$list_models[[1]]$nelson_allen_estimates)[2] - 1)
    #                                         for (i in 1:ntree){
    #                                           nelson_allen =
    #                                             nelson_allen +
    #                                             object$list_models[[i]]$nelson_allen_estimates[object$list_models[[i]]$nelson_allen_estimates[,1] == preds_node[j,i],
    #                                                                                2:ncol(object$list_models[[i]]$nelson_allen_estimates)]
    #                                         }
    #                                         return(nelson_allen / ntree)
    #                                       }))
    # Sys.time() - t1

    preds_surv_KMloc = exp(-mean_nelson_allen_estimates)

    predicted_KMloc =
      (preds_surv_KMloc[, which(object$list_models[[1]]$time < object$max_time)] -
         cbind(preds_surv_KMloc[, which(object$list_models[[1]]$time < object$max_time)[-1] ], 0)) %*%
      sapply(X = 1:length(c(object$list_models[[1]]$time[which(object$list_models[[1]]$time < object$max_time)[-1]], object$max_time)),
             FUN = function(i){do.call(object$phi,
                                       c(list(x=c(object$list_models[[1]]$time[which(object$list_models[[1]]$time < object$max_time)[-1]],
                                                  object$max_time)[i]),
                                         object$phi.args))})


    return(list(time = object$list_models[[1]]$time,
                preds_surv_KMloc = exp(-mean_nelson_allen_estimates),
                predicted_KMloc = as.vector(predicted_KMloc)))
  }
}




