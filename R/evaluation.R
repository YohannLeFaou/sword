
#' @title Compute the prediction of a model built with \code{\link{Cox_regression}}
#'
#' @description Given an output from \code{\link{Cox_regression}},
#' \code{predict_Cox_regression} allows to get its predictions for new
#' observations
#'
#' @param object A list output by \code{\link{Cox_regression}}
#' @param newdata A data.frame which contains the same variables as the ones
#' used for the training
#' @return A list with the following elements :
#' \item{predicted}{The vector of the predicted values for \code{phi}\eqn{(T')}
#' (with \eqn{T' = min(T, } \code{max_time}\eqn{)}) for the observations of \code{newdata}}
#' \item{survival}{The matrix which contains the estimated values of the survival curves at
#' \code{time_points}, for the observations of \code{newdata}}
#' \item{time_points}{The vector of the time points where the survival curves
#' are evaluated}
#'
#' @seealso \code{\link{Cox_regression}}
#'
#' @examples
#'
#' data(veteran, package = "randomForestSRC")
#' set.seed(17)
#' train_lines = sample(1:nrow(veteran), 120)
#' res1 = Cox_regression(y_var = "time",
#'                       delta_var = "status",
#'                       x_vars = setdiff(colnames(veteran),c("time","status")),
#'                       data_train = veteran[train_lines,],
#'                       types_weights_eval = c("KM", "Cox", "RSF", "unif"))
#'
#' pred1 = predict_Cox_regression(object = res1, newdata = veteran[-train_lines,])
#' print(pred1$predicted)


eval_aggregated_criteria = function(model_predictions, data, y_name, delta_name,
                                    max_time,  method, phi, phi.args,
                                    bandwidth = NULL){
  # this function computes the aggregated criteria for assessing performances of a model
  # parameters :
  ## model_predictions (vector) : predicted values of the model
  ## data : dataset with corresponding observations
  ## delta_name : name of the "delta" variable in the dataset "data"
  ## y_name : name of the "y" variable in the datasent "data"
  ## c_phi : threshold c_phi of the function phi

  data$pred = model_predictions
  data = data[order(data$pred, decreasing = F),]

  if (method == "group"){

    # add a column giving the group number
    data$group = (0:(length(model_predictions)-1) - (0:(length(model_predictions)-1) %% bandwidth))/ bandwidth + 1
    if(sum(data$group == max(data$group)) < bandwidth){
      data$group[which(data$group == max(data$group))] = max(data$group) - 1
    }
    n_groups = max(data$group)

    # compute empirical phi of the group & mean prediction inside the group
    phi_by_group = c()
    mean_pred_by_group = c()
    formula = stats::as.formula(paste0("Surv(time = ", y_name,", event = ", delta_name, ") ~ 1"))
    for (i in 1:n_groups){
      if (sum(data[data$group == i,delta_name]) > 0){ # if it exists non censored observations in the group

        # empirical phi inside the group (KM estimate)
        a = survival::survfit(formula = formula,
                              data = data[data$group == i,])
        phi_by_group[i] = - sum( diff(c(1,a$surv[which(a$time < max_time)],0)) *
                                   sapply(X = 1:length(c(a$time[which(a$time < max_time)], max_time)),
                                          FUN = function(i){do.call(phi, c(list(x=c(a$time[which(a$time < max_time)], max_time)[i]),
                                                                           phi.args))})
        )

        # mean prediction for the group
        mean_pred_by_group[i] = mean(data$pred[data$group == i])

      } else { # if all observations in the group are censored
        phi_by_group[i] = NA
        mean_pred_by_group[i] = NA
      }
    }
    phi_by_group = phi_by_group[!is.na(phi_by_group)]
    mean_pred_by_group = mean_pred_by_group[!is.na(mean_pred_by_group)]

    # compare empirical phi of the group with mean prediction inside the group
    return(c(NormalizedGini(solutions = 1:length(phi_by_group), # Gini
                            predictions = phi_by_group),
             (1 + stats::cor.test(1:length(phi_by_group), phi_by_group, method = "kendall")$estimate)/2, # Kendall
             (1 - 1/length(phi_by_group) * sum((phi_by_group - mean_pred_by_group)^2) / stats::sd(phi_by_group)^2) # R2
    ))
  }
}

#' @title helper function
#'
#' @description function called by \code{\link{NormalizedGini}}
#'
#' @seealso \code{\link{NormalizedGini}}


SumModelGini <- function(solutions, predictions, weights){
  ## function called by NormalizedGini

  df = data.frame(rank_solutions = rank(solutions), predictions = predictions, weights = weights)
  df <- df[order(df$predictions, decreasing = TRUE),]
  df$random = (1:nrow(df))/nrow(df)
  totalPos <- sum(df$weights * df$rank_solutions)
  df$cumPosFound <- cumsum(df$weights * df$rank_solutions) # this will store the cumulative number of positive examples found (used for computing "Model Lorentz")
  df$Lorentz <- df$cumPosFound / totalPos # this will store the cumulative proportion of positive examples found ("Model Lorentz")
  df$Gini <- df$Lorentz - df$random # will store Lorentz minus random
  return(sum(df$Gini))
}


#' @title Compute a Gini goodness of fit statistic
#'
#' @description Given a vector of observed values \code{solutions} and a vector
#' of predicted values \code{predictions}, the Gini index
#' measures how well the order of the predicted values corresponds
#' to the order of the observed values.
#'
#' @param solutions A vector of observed values
#' @param predictions A vector of predicted values
#' @param weights A vector of weights for the single observations (défault = \code{NULL}).
#' If \code{NULL}, then weights are taken as equal to 1
#'
#' @details
#' \itemize{
#' \item \code{solutions} is transformed into \code{rank(solutions)} before we compute
#' the gini coefficient
#' \item The weighted Gini index may not give meaningfull information if \code{weights}
#' depends on \code{predictions}
#' }
#'
#' @return The real number giving the value of the Gini index
#'
#' @seealso \code{\link{SumModelGini}}
#'
#' @examples
#'


NormalizedGini <- function(solutions, predictions, weights = NULL) {
  # function which computes the Gini index of performance of a model
  # solutions : vector of true values
  # predictions : vector of predicted values

  if(is.null(weights)){weights = rep(1, length(solutions))}
  SumModelGini(solutions, predictions, weights) / SumModelGini(solutions, solutions, weights)
}


#' @title Compute weighted quadratic errors statistics
#'
#' @description Given a vector of observed values \code{solutions} and a vector
#' of predicted values \code{predictions}, compute the Mean Squared Error (MSE) and
#' the percentage of explained variance (R2) statistics.
#'
#' @param solutions A vector of observed values
#' @param predictions A vector of predicted values
#' @param weights A vector of weights for the single observations (défault = \code{NULL}).
#' If \code{NULL}, then weights are taken as equal to 1
#'
#' @return The vector with the value of the MSE (index 1) and the R2 statistics
#'
#' @seealso \code{\link{SumModelGini}}
#'
#' @examples
#'

eval_weighted_criteria = function(predictions, solutions, weights = NULL){
  if(!is.null(weights)){weights = weights / sum(weights)}
  if(is.null(weights)){weights = rep(1 / length(solutions), length(solutions))}
  mean = sum(weights * solutions)
  scores = c(sum(weights * (solutions - predictions)^2), # mse
             (1 - sum(weights * (solutions - predictions)^2) / # R2
                sum(weights * (solutions - mean)^2))
  )
  names(scores) = c("mse", "R2")
  return(scores)
}


eval_model = function(predictions,
                      data,
                      phi_name,
                      y_name,
                      delta_name,
                      max_time,
                      eval_methods,
                      phi = phi,
                      phi.args,
                      mat_weights = NULL,
                      phi_non_censored_name = NULL,
                      v_bandwidth){ # pour le no_w je peux mettre un vecteur avec 0 et 1 dans no_w

  # rmq : à priori on pourrait se passer d'avoir "data" dans les arguments
  # (et se contenter de phi, y(_prime), delta(_prime))
  # mais je laisse au cas où

  list_criteria = list()

  if("concordance" %in% eval_methods){
    formula = stats::as.formula(paste0("Surv(", phi_name, ",", delta_name, ") ~ predictions"))
    concordance =
      1 - survival::survConcordance(formula = formula,
                                    data = cbind(data[,c(phi_name, delta_name)],
                                                 predictions)
      )$concordance
    list_criteria$concordance = concordance
  }

  if ("group" %in% eval_methods){
    criteria_group = as.vector(sapply(X = v_bandwidth,
                                      FUN = eval_aggregated_criteria,
                                      model_predictions = predictions,
                                      data = data[,c("y_prime","delta_prime")],
                                      y_name = "y_prime",
                                      delta_name = "delta_prime",
                                      max_time = max_time,
                                      method = "group",
                                      phi = phi,
                                      phi.args = phi.args
    ))
    names(criteria_group) = c(do.call(c, lapply( X = v_bandwidth, FUN = function(x){paste0(c("Gini_","Kendall_","R2_"), x)})))
    list_criteria$criteria_group = criteria_group
  }

  if ("weighted" %in% eval_methods){
    criteria_weighted = as.vector(apply(X = mat_weights,
                                        MARGIN = 2,
                                        FUN = eval_weighted_criteria,
                                        predictions = predictions,
                                        solutions = data[,phi_name]))

    names(criteria_weighted) = as.vector(sapply(X = colnames(mat_weights),
                                                FUN = function(x){paste0(x, c("_mse", "_R2"))}))
    list_criteria$criteria_weighted = criteria_weighted
  }

  if (!is.null(phi_non_censored_name)){

    mse = 1/length(data[,phi_non_censored_name]) *
      sum( (data[,phi_non_censored_name] - predictions)^2 )

    R2 = 1 - sum( (data[,phi_non_censored_name] - predictions)^2 ) /
      sum( (data[,phi_non_censored_name] - mean(data[,phi_non_censored_name]))^2 )

    Kendall = (1 + stats::cor.test(predictions,
                            data[,phi_non_censored_name],
                            method = "kendall")$estimate)/2

    Gini = NormalizedGini(predictions,
                          data[,phi_non_censored_name])

    criteria_non_censored = c(Gini, Kendall, R2, mse)
    names(criteria_non_censored) = c("Gini", "Kendall", "R2", "mse")
    list_criteria$criteria_non_censored = criteria_non_censored
  }
  return(list_criteria)
}


