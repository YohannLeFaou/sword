

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
    formula = as.formula(paste0("survival::Surv(time = ", y_name,", event = ", delta_name, ") ~ 1"))
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
             (1 + cor.test(1:length(phi_by_group), phi_by_group, method = "kendall")$estimate)/2, # Kendall
             (1 - 1/length(phi_by_group) * sum((phi_by_group - mean_pred_by_group)^2) / sd(phi_by_group)^2) # R2
    ))
  }
  if (method == "single"){

    # add a column giving the group number
    index_event = which(data[,delta_name] == 1)
    n_groups = length(index_event)
    if (data[1,delta_name] == 0){
      groups = rep(x = 1:n_groups, times = index_event - c(0,index_event[1:(n_groups-1)]))
      data = data[1:length(groups),]
    }
    if (data[1,delta_name] == 1){
      groups = rep(x = 1:(n_groups-1), times = (index_event[2:(n_groups)] ) - c(index_event[1:(n_groups-1)]))
      groups = c(groups, rep(n_groups, times = nrow(data) - length(groups) ))
    }
    data$group = groups

    # compute empirical phi of the group & mean prediction inside the group
    phi_by_group = c()
    mean_pred_by_group = c()
    weights = c()
    for (i in 1:n_groups){

      # empirical phi inside the group (KM estimate)
      b = data[which(data$group == i),]
      jump_size = 1 / sum(b[,y_name] >= b[which(b[,delta_name] == 1), y_name])
      phi_by_group[i] =
        jump_size * do.call(phi, c(list(x=b[which(b[,delta_name] == 1),y_name]), phi.args)) +
        (1 - jump_size) * do.call(phi, c(list(x=max_time), phi.args))

      mean_pred_by_group[i] = mean(data$pred[data$group == i])
      weights[i] = nrow(b)
    }

    # compare empirical phi of the group with mean prediction inside the group
    mean_phi_by_group = sum( weights * phi_by_group ) / sum(weights)

    return(c(NormalizedGini(solutions = 1:sum(weights),
                            predictions = rep(x = phi_by_group, times = weights)),
             (1 + cor.test(1:sum(weights),
                           rep(x = phi_by_group, times = weights),
                           method = "kendall")$estimate)/2,
             (1 - sum( weights * (phi_by_group - mean_pred_by_group)^2) /
                sum( weights * (phi_by_group - mean_phi_by_group)^2))
    ))
  }
}


SumModelGini <- function(solutions, predictions, weights){
  ## function called by NormalizedGini

  df = data.frame(solutions = solutions, predictions = predictions, weights = weights)
  df <- df[order(df$predictions, decreasing = TRUE),]
  df$random = (1:nrow(df))/nrow(df)
  totalPos <- sum(df$weights * df$solutions)
  df$cumPosFound <- cumsum(df$weights * df$solutions) # this will store the cumulative number of positive examples found (used for computing "Model Lorentz")
  df$Lorentz <- df$cumPosFound / totalPos # this will store the cumulative proportion of positive examples found ("Model Lorentz")
  df$Gini <- df$Lorentz - df$random # will store Lorentz minus random
  return(sum(df$Gini))
}

NormalizedGini <- function(solutions, predictions, weights = NULL) {
  # function which computes the Gini index of performance of a model
  # solutions : vector of true values
  # predictions : vector of predicted values

  if(is.null(weights)){weights = rep(1, length(solutions))}
  SumModelGini(solutions, predictions, weights) / SumModelGini(solutions, solutions, weights)
}



eval_weighted_criteria = function(predictions, solutions, weights){
  mean = sum(weights * solutions)
  return(
    c(sum(weights * (solutions - predictions)^2), # mse
      (1 - sum(weights * (solutions - predictions)^2) / # R2
         sum(weights * (solutions - mean)^2))
    )
  )
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
    formula = as.formula(paste0("survival::Surv(", phi_name, ",", delta_name, ") ~ predictions"))
    concordance =
      1 - survival::survConcordance(formula = formula,
                                    data = cbind(data[,c(phi_name, delta_name)],
                                                 predictions)
      )$concordance
    list_criteria$concordance = concordance
  }


  if ("single" %in% eval_methods){
    criteria_single = eval_aggregated_criteria(model_predictions = predictions,
                                               data = data[,c("y_prime","delta_prime")],
                                               y_name = "y_prime",
                                               delta_name = "delta_prime", # important to set delta_prime such the group are good in method "single"
                                               max_time = max_time,
                                               method = "single",
                                               phi = phi,
                                               phi.args = phi.args
    )
    names(criteria_single) = c("Gini_single","Kendall_single","R2_single")
    list_criteria$criteria_single = criteria_single
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

    Kendall = (1 + cor.test(predictions,
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


