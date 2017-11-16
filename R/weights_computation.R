

make_weights_from_survival_curves = function(vect_y, vect_delta, mat_survival_curves_C, time_points){
  # function which computes weights of the observations given values of y, delta and estimates
  ## of S_C(.|X)

  x = cbind(vect_y, mat_survival_curves_C)
  x = x[(vect_delta == 1),]
  non_zero_w = 1 / length(vect_delta) /
    apply(X = x,
          MARGIN = 1,
          FUN = function(a){stats::approx(x = time_points,
                                   y = a[-1],
                                   xout = a[1],
                                   rule = 2,
                                   method = "linear")$y
          }
    )
  w = vect_delta
  w[w != 0] = non_zero_w
  return(w)
}


make_KM_weights = function(vect_y, vect_delta, vect_c = NULL){
  # function which computes the weights of the observations in the Kaplan Meier case
  ## vect_c : values taken by C in the case C is never censored

  is_C_censored = is.null(vect_c)
  # is_C_censored = FALSE if C is not censored (in some situations, we allways observe C)
  # is_C_censored = TRUE if C is censored (by the time of interest)

  if(is_C_censored){
    # S_C is estimated by Kaplan Meier estimator
    delta_c = 1 - vect_delta
    data_c = data.frame(y = vect_y, delta_c = delta_c)
    km_c = survival::survfit(survival::Surv(time = y, event = delta_c) ~ 1,
                             data = data_c,
                             type = "kaplan-meier")
  } else{
    # S_C is estimated by the empirical estimator of CDF
    delta_c = rep(1, length(vect_y))
    data_c = data.frame(y = vect_c,
                        delta_c = delta_c)
    km_c = survival::survfit(survival::Surv(time = y , event = delta_c) ~ 1,
                             data = data_c,
                             type = "kaplan-meier")
  }
  c_time = c(0,km_c$time)
  c_surv = c(1,km_c$surv)
  N = length(vect_y)
  w = stats::approx(x = c_time, # weights are estimated by interpolation of S_C estimator
             y = c_surv,
             xout = vect_y,
             method = "linear",
             rule = 2)
  weights = ifelse(vect_delta == 0, 0, 1/N * 1/w$y)
  return(weights)
}


make_weights = function(data,
                        y_name,
                        delta_name,
                        y_name2 = NULL,
                        delta_name2 = NULL,
                        type = c("KM", "COX", "RSF", "0_1"),
                        max_ratio_weights,
                        x_vars = NULL,
                        censoring_model_object = TRUE){

  type = match.arg(as.character(type), c("KM", "Cox", "RSF", "0_1"))
  if (is.null(y_name2)) y_name2 = y_name
  if (is.null(delta_name2)) delta_name2 = delta_name

  data$deltaC = 1 * (data[,delta_name2] == 0)
  weight_max = max_ratio_weights / nrow(data)
  list_result = list()

  # rmk : no need to fill list_result$censoring_model_object if "censoring_model_object = TRUE" :
  # list_result$censoring_model_object will return "NULL" by default

  if (type == "0_1"){
    weights = data[,delta_name]
  }

  if (type == "KM"){
    if (sum(data[,delta_name]) == nrow(data)){
      weights = rep(1/nrow(data), nrow(data))
    } else {
      weights = make_KM_weights(vect_y = data[,y_name],
                                vect_delta = data[,delta_name])
    }
  }

  if (type == "Cox"){
    # possible to make a functions "make_cox_weights", "make_RSF_weights"
    formula = stats::as.formula(paste0("survival::Surv(",y_name2,", deltaC ) ~ ."))
    cox_fit = survival::coxph(formula = formula,
                              data = data[,c(y_name2, "deltaC", x_vars)]
    )

    baseline_cox = survival::basehaz(cox_fit)
    ref_surv = data.frame(time = c(0, baseline_cox$time),
                          surv = c(1, exp(-baseline_cox$hazard)))

    weights = 1/ nrow(data) /
      pmax( (stats::approx(x = ref_surv$time,
                    y = ref_surv$surv,
                    xout = data[,y_name],
                    method = "linear",
                    rule = 2)$y)^(exp(cox_fit$linear.predictors)) , 1/(max_ratio_weights * 1.5) )
    weights[data[,delta_name] == 0] = 0
    if(censoring_model_object){list_result$censoring_model_object = cox_fit}
  }

  if (type == "RSF"){
    max_time = max(data[which(data[,delta_name] == 1), y_name])
    ntime = seq(from = 0, to = max_time * 1.05, length.out = 100)
    nodedepth_RSF = floor(1/2*log(nrow(data))/log(2))
    formula = stats::as.formula(paste0("survival::Surv(",y_name2,", deltaC ) ~ ."))

    RSF_fit = randomForestSRC::rfsrc(formula = formula,
                                     data = data[, c(y_name2, "deltaC", x_vars)],
                                     forest = T,
                                     ntree = 100,
                                     nsplit = 10,
                                     nodedepth = nodedepth_RSF,
                                     splitrule = "logrank",
                                     ntime = ntime)

    weights = make_weights_from_survival_curves(vect_y = data[, y_name],
                                                vect_delta = data[, delta_name],
                                                mat_survival_curves_C = cbind(1, RSF_fit$survival.oob),
                                                time_points = c(0, RSF_fit$time.interest))
    if(censoring_model_object){list_result$censoring_model_object = RSF_fit}
  }
  sum_weights = sum(weights)
  n_weights_modif = sum(weights > (min(weights[weights > 0]) *  max_ratio_weights))
  weights = pmin(weights, (min(weights[weights > 0]) *  max_ratio_weights) )

  list_result$weights = weights / sum(weights)
  list_result$sum_weights = sum_weights
  list_result$n_weights_modif = n_weights_modif

  return(list_result)
}

