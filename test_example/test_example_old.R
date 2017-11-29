
# ------------------------------------------
#           Cox_regression
# ------------------------------------------

# ------------------------------------------------
#   Basic call to train a model
# ------------------------------------------------

data(veteran, package = "randomForestSRC")
res1 = Cox_regression(y_var = "time",
                      delta_var = "status",
                      x_vars = setdiff(colnames(veteran),c("time","status")),
                      data_train = veteran,
                      types_weights_eval = c("KM", "Cox", "RSF", "0_1"))

print(res1$list_criteria_train)
matplot(y = t(res1$survival_train), x = res1$time_points, type = "l")

print(res1$max_time) # by default \code{max_time} is set to 999 which is too large
# to have good predictions (train error underestimate the error of the model)

# ------------------------------------------------
#   Training with estimation of test error
# ------------------------------------------------
set.seed(17)
train_lines = sample(1:nrow(veteran), 100)
res2 = Cox_regression(y_var = "time",
                      delta_var = "status",
                      x_vars = setdiff(colnames(veteran),c("time","status")),
                      data_train = veteran[train_lines,],
                      data_test = veteran[-train_lines,],
                      types_weights_eval = c("KM", "Cox", "RSF", "0_1"))

print(res2$max_time)
print(res2$list_criteria_train,5)
print(res2$list_criteria_test, 5) # test error is low
print(res2$predicted_test)
matplot(y = t(res2$survival_test), x = res2$time_points, type = "l")

# ------------------------------------------------
#   Modify the max_time
# ------------------------------------------------
set.seed(17)
train_lines = sample(1:nrow(veteran), 100)
res3 = Cox_regression(y_var = "time",
                      delta_var = "status",
                      x_vars = setdiff(colnames(veteran),c("time","status")),
                      data_train = veteran[train_lines,],
                      data_test = veteran[-train_lines,],
                      max_time = 600,
                      types_weights_eval = c("KM", "Cox", "RSF", "0_1"))

print(res3$list_criteria_train)
print(res3$list_criteria_test) # test error is much better
print(res3$predicted_test)
matplot(y = t(res3$survival_test), x = res3$time_points, type = "l")

# ------------------------------------------------
#   Use custom \code{phi} function
# ------------------------------------------------
g = function(x,a) abs(x-a)
set.seed(17)
train_lines = sample(1:nrow(veteran), 100)
res4 = Cox_regression(y_var = "time",
                      delta_var = "status",
                      x_vars = setdiff(colnames(veteran),c("time","status")),
                      data_train = veteran[train_lines,],
                      data_test = veteran[-train_lines,],
                      phi = g,
                      phi.args = list(a = 200),
                      max_time = 400,
                      types_weights_eval = c("KM", "Cox", "RSF", "0_1")
                      )

print(res4$list_criteria_test)
print(res4$predicted_test)
print(res4$censoring_rate_with_threshold) # very low rate of censoring
print(head(res4$mat_weights_test))
print(apply(X = res4$mat_weights_test, MARGIN = 2, FUN = function(x){max(x[x != 0])/min(x[x != 0])}))
# ------------------------------------------------
#   Modify \code{max_ratio_weights_eval}
# ------------------------------------------------

set.seed(17)
train_lines = sample(1:nrow(veteran), 100)
res5 = Cox_regression(y_var = "time",
                      delta_var = "status",
                      x_vars = setdiff(colnames(veteran),c("time","status")),
                      data_train = veteran[train_lines,],
                      data_test = veteran[-train_lines,],
                      phi = g,
                      phi.args = list(a = 200),
                      max_time = 400,
                      max_ratio_weights_eval = 5,
                      types_weights_eval = c("KM", "Cox", "RSF", "0_1")
)

print(res4$list_criteria_test)
print(res4$predicted_test)

# ------------------------------------------
#           predict_Cox_regression
# ------------------------------------------
data(veteran, package = "randomForestSRC")
set.seed(17)
train_lines = sample(1:nrow(veteran), 120)
res1 = Cox_regression(y_var = "time",
                      delta_var = "status",
                      x_vars = setdiff(colnames(veteran),c("time","status")),
                      data_train = veteran[train_lines,],
                      types_weights_eval = c("KM", "Cox", "RSF", "0_1"))

pred1 = predict_Cox_regression(object = res1, newdata = veteran[-train_lines,])
print(pred1$predicted)


# ------------------------------------------
#           RSF_regression
# ------------------------------------------
data(veteran, package = "randomForestSRC")
res1 = RSF_regression(y_var = "time",
                      delta_var = "status",
                      x_vars = setdiff(colnames(veteran),c("time","status")),
                      data_train = veteran,
                      types_weights_eval = c("KM", "Cox", "RSF", "0_1"))
print(res1$list_criteria_train)

set.seed(17)
train_lines = sample(1:nrow(veteran), 100)
res2 = RSF_regression(y_var = "time",
                      delta_var = "status",
                      x_vars = setdiff(colnames(veteran),c("time","status")),
                      data_train = veteran[train_lines,],
                      data_test = veteran[-train_lines,],
                      phi = sqrt,
                      max_time = 300,
                      types_weights_eval = c("KM", "Cox", "RSF", "0_1"),
                      nodedepth = 3 #additionnal parameter provided for the random survival
                                    #forest (randomForestSRC)
                      )
print(res2$list_criteria_train)
print(res2$list_criteria_test)
print(res2$predicted_test)
matplot(y = t(res2$survival_test), x = res2$time_points, type = "l")


# --------------------------------------------
#           weighted_regression_survival
# --------------------------------------------

set.seed(17)
train_lines = sample(1:nrow(veteran), 100)
res2 = RSF_regression(y_var = "time",
                      delta_var = "status",
                      x_vars = setdiff(colnames(veteran),c("time","status")),
                      data_train = veteran[train_lines,],
                      data_test = veteran[-train_lines,],
                      phi = sqrt,
                      max_time = 300,
                      types_weights_eval = c("KM", "Cox", "RSF", "0_1"),
                      nodedepth = 3 #additionnal parameter provided for the random survival
                      #forest (randomForestSRC)
)
print(res2$list_criteria_train)
print(res2$list_criteria_test)
print(res2$predicted_test)
matplot(y = t(res2$survival_test), x = res2$time_points, type = "l")

weighted_regression_survival()


####################################################
####### autres tests à supprimer dans le futur

make_KM_weights(vect_y = veteran$time,
                vect_delta = veteran$status)

make_weights(data = veteran, y_name = "time", delta_name = "status", type = "Cox", max_ratio_weights = 1000,
             x_vars = setdiff(colnames(veteran),c("time","status")), censoring_model_object = T)

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

predict_nodes2 = function (object, newdata, na.action = stats::na.pass) {
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
    rpart:::pred.rpart(object, as.matrix(newdata))
    predict(object, newdata, type = "vector")
  }
  as.integer(row.names(object$frame))[where]
}


library(rpart)
fit <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
predict_nodes(object = fit, newdata = kyphosis, na.action = stats::na.pass)
predict_nodes2(object = fit, newdata = kyphosis, na.action = stats::na.pass)

