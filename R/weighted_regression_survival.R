

weighted_regression_survival = function(time_var,
                                        event_var,
                                        x_vars,
                                        data_train,
                                        data_test = NULL,
                                        type_regression = c("RF", "gam"),
                                        type_weights = c("KM","Cox","RSF"),
                                        phi = function(x){x},
                                        phi.args = list(),
                                        max_time = NULL,
                                        weighted_regression_object = T,
                                        censoring_model_object = T,
                                        eval_methods = c("concordance","weighted"),
                                        v_bandwidth = 20,
                                        types_weights_eval = c("KM"),
                                        max_ratio_weights_model = NULL,
                                        max_ratio_weights_eval = 20,
                                        mat_weights = NULL, # we say that when mat_weights is not NULL, and type_weights is NULL, the column 1 of mat_weights is employed for fitting
                                        time_non_censored_var = NULL,
                                        mode_w_RF = 1,
                                        ...){

  # Preprocessing of the arguments & data

  ## specific preprocessing compared to Cox/RSF regression
  type_weights = match.arg(as.character(type_weights), c("KM","Cox","RSF"))
  weights_manual = !is.null(mat_weights)
  type_regression = match.arg(as.character(type_regression), c("RF", "gam"))
  # --------------------------------------------------------------------------------

  eval_methods <- match.arg(as.character(eval_methods), c("concordance","single", "group", "weighted"), several.ok = T)
  types_weights_eval = match.arg(as.character(types_weights_eval), c("KM", "Cox", "RSF", "0_1"), several.ok = T)
  types_weights_eval = unique(c(type_weights, types_weights_eval)) # weights for trainig are used for evaluation

  # column names of mat_weights should be explicit
  if (!is.null(mat_weights) & is.null(colnames(mat_weights))) colnames(mat_weights) = paste0("w",1:ncol(mat_weights))

  if (is.null(time_non_censored_var)) {
    phi_non_censored_name = NULL
  } else {
    phi_non_censored_name = "phi_non_censored"
  }

  if (!is.null(data_test)){
    data = rbind(data_train[,c(time_var, event_var, x_vars, time_non_censored_var)],
                 data_test[,c(time_var, event_var, x_vars, time_non_censored_var)])
    data$is_train = c(rep(1, nrow(data_train)), rep(0, nrow(data_test)))
  } else {
    data = data_train[,c(time_var, event_var, x_vars, time_non_censored_var)]
    data$is_train = 1
  }

  if (is.null(max_time)){max_time = max(data_train[which(data_train[, event_var] == 1), time_var])}

  data$y_prime = pmin(data[,time_var], max_time)
  data$delta_prime = 1 * ((data[,event_var] != 0) | (data[,time_var] >= max_time))
  data$phi = sapply(X = 1:length(data$y_prime),
                    FUN = function(i){do.call(phi, c(list(x=data$y_prime[i]), phi.args))})
  if(!is.null(time_non_censored_var)){
    data$phi_non_censored = sapply(X = 1:nrow(data),
                                   FUN = function(i){do.call(phi, c(list(x=pmin(data[,time_non_censored_var], max_time)[i]), phi.args))})
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
      v_weights_model_test = mat_weights[(nrow(data_train) +1):nrow(data), type_weights]
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

    for (j in 1:length(types_weights_eval)){
      if (types_weights_eval[j] == type_weights){
        # need to differentiate type_weights or not because of the option "censoring_model_object"
        res_weights_train = make_weights(data = data[data$is_train == 1, ],
                                         y_name = "y_prime",
                                         delta_name = "delta_prime",
                                         y_name2 = time_var,
                                         delta_name2 = event_var,
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
                                          y_name2 = time_var,
                                          delta_name2 = event_var,
                                          type = types_weights_eval[j],
                                          max_ratio_weights = 1000,
                                          x_vars = x_vars,
                                          censoring_model_object = FALSE)
          mat_weights_test[,j] = res_weights_test$weights
          v_weights_model_test = res_weights_test$weights
        }
      } else {
        res_weights_train = make_weights(data = data[data$is_train == 1, ],
                                         y_name = "y_prime",
                                         delta_name = "delta_prime",
                                         y_name2 = time_var,
                                         delta_name2 = event_var,
                                         type = types_weights_eval[j],
                                         max_ratio_weights = 1000,
                                         x_vars = x_vars,
                                         censoring_model_object = FALSE)
        mat_weights_train[,j] = res_weights_train$weights

        if (!is.null(data_test)){
          res_weights_test = make_weights(data = data[data$is_train == 0, ],
                                          y_name = "y_prime",
                                          delta_name = "delta_prime",
                                          y_name2 = time_var,
                                          delta_name2 = event_var,
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

  ## test
  if (!is.null(data_test)){
    n_weights_model_modif_test = sum(v_weights_model_test > (min(v_weights_model_test[v_weights_model_test > 0]) * max_ratio_weights_model))
    v_weights_model_test = pmin(v_weights_model_test, min(v_weights_model_test[v_weights_model_test > 0]) * max_ratio_weights_model)
    v_weights_model_test = v_weights_model_test / sum(v_weights_model_test)
  }

  # weights_eval
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

  # build train & test
  data_train = data[data$is_train == 1,]
  if (!is.null(data_test)){
    data_test = data[data$is_train == 0,]
  }


  weighted_regression_result = make_res_weighted_regression(time_var = time_var,
                                                            event_var = event_var,
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
                                                            ...)

  weighted_regression_result$type_weights = type_weights
  weighted_regression_result$max_ratio_weights_model = max_ratio_weights_model
  weighted_regression_result$mode_w_RF = mode_w_RF
  weighted_regression_result$sum_weights_train = sum_weights_train
  weighted_regression_result$n_weights_model_modif_train = n_weights_model_modif_train
  if (!is.null(data_test)){
    weighted_regression_result$sum_weights_test = sum_weights_test
    weighted_regression_result$n_weights_model_modif_test = n_weights_model_modif_test
  }
  if (!is.null(censoring_model_object)){
    # return only the model for censoring fitted on train data
    weighted_regression_result$censoring_model_object = censoring_model_object
  }
  return(weighted_regression_result)
}


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
    predictions = predict_rpartRF(object = object$weighted_rpartRF_object,
                                  newdata = newdata[,object$x_vars],
                                  type = "response")
    predictions_surv = predict_rpartRF(object = object$weighted_rpartRF_object,
                                       newdata = newdata[,object$x_vars],
                                       type = "surv")
    return(list(predicted = predictions,
                predicted_surv = predictions_surv))
  }
}


make_res_weighted_regression = function(time_var,
                                        event_var,
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
                                        ...){

  # Build and predict on train
  if (type_regression == "RF"){
    if (mode_w_RF == 1){
      RF_fit = randomForestSRC::rfsrc(formula = phi ~ .,
                                      data = data_train[v_weights_model_train > 0, c("phi",x_vars)],
                                      case.wt =  v_weights_model_train[v_weights_model_train > 0],
                                      forest = T,
                                      ...)

      overfitted_predictions = randomForestSRC::predict.rfsrc(RF_fit,
                                                              data_train[,x_vars])$predicted
    }
    if (mode_w_RF == 2){
      rpartRF_fit = rpartRF(data = data_train,
                            time_var = time_var,
                            event_var = event_var,
                            x_vars = x_vars,
                            type_weights = type_weights,
                            max_ratio_weights_model = max_ratio_weights_model,
                            ...)

      overfitted_predictions = predict_rpartRF(object = rpartRF_fit,
                                               newdata = data_train[,x_vars],
                                               type = "response")


      res_overfitted_surv_curv_KMloc = predict_rpartRF(object = rpartRF_fit,
                                                       newdata = data_train[,x_vars],
                                                       type = "surv")


      overfitted_predictions_KMloc =
        (res_overfitted_surv_curv_KMloc$surv[, which(res_overfitted_surv_curv_KMloc$time < max_time)] -
           cbind(res_overfitted_surv_curv_KMloc$surv[, which(res_overfitted_surv_curv_KMloc$time < max_time)[-1] ], 0)) %*%
        sapply(X = 1:length(c(res_overfitted_surv_curv_KMloc$time[which(res_overfitted_surv_curv_KMloc$time < max_time)[-1]], max_time)),
               FUN = function(i){do.call(phi,
                                         c(list(x=c(res_overfitted_surv_curv_KMloc$time[which(res_overfitted_surv_curv_KMloc$time < max_time)[-1]],
                                                    max_time)[i]),
                                           phi.args))})

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
                                           type = "response")

        res_test_surv_curv_KMloc = predict_rpartRF(object = rpartRF_fit,
                                                   newdata = data_test[,x_vars],
                                                   type = "surv")

        test_predictions_KMloc =
          (res_test_surv_curv_KMloc$surv[, which(res_test_surv_curv_KMloc$time < max_time)] -
             cbind(res_test_surv_curv_KMloc$surv[, which(res_test_surv_curv_KMloc$time < max_time)[-1] ], 0)) %*%
          sapply(X = 1:length(c(res_test_surv_curv_KMloc$time[which(res_test_surv_curv_KMloc$time < max_time)[-1]], max_time)),
                 FUN = function(i){do.call(phi,
                                           c(list(x=c(res_test_surv_curv_KMloc$time[which(res_test_surv_curv_KMloc$time < max_time)[-1]],
                                                      max_time)[i]),
                                             phi.args))})
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
    data_train = data_train[,c(time_var, event_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)],
    v_weights_model_train = v_weights_model_train,
    mat_weights_train = mat_weights_train,
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
    result$res_surv_curv_train_KMloc = res_overfitted_surv_curv_KMloc
    result$predicted_train_KMloc = overfitted_predictions_KMloc
  }
  if (!is.null(data_test)){
    result$predicted_test = test_predictions
    result$list_criteria_test = list_criteria_test
    result$data_test = data_test[,c(time_var, event_var, "y_prime", "delta_prime", "phi", phi_non_censored_name, x_vars)]
    result$mat_weights_test = mat_weights_test
    if ((type_regression == "RF") & (mode_w_RF == 2)){
      result$list_criteria_test_KMloc = list_criteria_test_KMloc
      result$res_surv_curv_test_KMloc = res_test_surv_curv_KMloc
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
                   time_var,
                   event_var,
                   x_vars,
                   type_weights,
                   max_ratio_weights_model,
                   ntree, # peut-etre a enlever plus tard
                   ...){
  list_models = list()
  for (i in 1:ntree){
    cat(paste0(i," "))
    sample_train = sample(x = 1:nrow(data), size = nrow(data), replace = T)
    d_train = data[sample_train, ]

    weights_in_bag = make_weights(y_name = "y_prime",
                                  delta_name = "delta_prime",
                                  y_name2 = time_var,
                                  delta_name2 = event_var,
                                  x_vars = x_vars,
                                  censoring_model_object = F,
                                  max_ratio_weights = max_ratio_weights_model,
                                  type = type_weights,
                                  data = d_train)$weights

    model = rpart::rpart(formula = phi ~ .,
                         data = d_train[which(weights_in_bag > 0), c("phi", x_vars)],
                         weights = weights_in_bag[which(weights_in_bag > 0)],
                         ...)

    pred_nodes = predict_nodes(object = model, newdata = d_train[,x_vars])
    nelson_allen_estimates = do.call(rbind,
                                     args = lapply(X = unique(pred_nodes),
                                                   FUN = function(node){
                                                     fit = survival::survfit(formula = stats::as.formula(paste0("Surv(time = ", time_var,", event = ",  event_var, ") ~ 1")),
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
  return(list_models)
}

predict_rpartRF = function(object, newdata, type){
  if (type == "response"){
    preds = do.call(what = cbind, args = lapply(X = object, FUN = function(x){rpart:::predict.rpart(object = x$model, newdata = newdata)}))
    return(rowMeans(preds))
  }
  if (type == "surv"){
    preds_node = do.call(what = cbind, args = lapply(X = object, FUN = function(x){predict_nodes(object = x$model, newdata = newdata)}))
    ntree = ncol(preds_node)

    #not as fast as the next solution but no bug here
    #t1 = Sys.time()
    mean_nelson_allen_estimates = do.call(what = rbind,
            args = lapply(X = 1:dim(preds_node)[1], FUN = function(j){ # t(preds_node)
              mat_nelson_allen = do.call(what = rbind,
                      args = lapply(X = 1:ntree, FUN = function(i){
                        object[[i]]$nelson_allen_estimates[object[[i]]$nelson_allen_estimates[,1] == preds_node[j,i],
                                                           2:ncol(object[[i]]$nelson_allen_estimates)]
                      }))
              colMeans(mat_nelson_allen)
    }))
    # #Sys.time() - t1

    # t1 = Sys.time()
    # mean_nelson_allen_estimates = do.call(what = rbind,
    #                                       args = lapply(X = 1:dim(preds_node)[1], FUN = function(j){
    #                                         nelson_allen = rep(0, dim(object[[1]]$nelson_allen_estimates)[2] - 1)
    #                                         for (i in 1:ntree){
    #                                           nelson_allen =
    #                                             nelson_allen +
    #                                             object[[i]]$nelson_allen_estimates[object[[i]]$nelson_allen_estimates[,1] == preds_node[j,i],
    #                                                                                2:ncol(object[[i]]$nelson_allen_estimates)]
    #                                         }
    #                                         return(nelson_allen / ntree)
    #                                       }))
    # Sys.time() - t1

    return(list(time = object[[1]]$time, surv = exp(-mean_nelson_allen_estimates)))
  }
}

