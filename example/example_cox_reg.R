
# ------------------------------------------------------------------------------------
#
#                             cox_reg
#
# ------------------------------------------------------------------------------------

# ------------------------------------------------
#   Load "transplant" data
# ------------------------------------------------
data("transplant", package = "survival")
transplant$delta = 1 * (transplant$event == "ltx") # create binary var
# which indicate censoring/non censoring

# keep only rows with no missing value
apply(transplant, MARGIN = 2, FUN = function(x){sum(is.na(x))})
transplant_bis = transplant[stats::complete.cases(transplant),]

# plot the survival curve of transplant data
KM_transplant = survfit(formula = survival::Surv(time = futime, event = delta) ~ 1,
                        data = transplant_bis)
plot(KM_transplant)

# ------------------------------------------------
#   Basic call to train a model
# ------------------------------------------------

res1 = cox_reg(y_var = "futime",
               delta_var = "delta",
               x_vars = setdiff(colnames(transplant_bis),
                                c("futime", "delta", "event")),
               train = transplant_bis,
               types_w_ev = c("KM", "Cox", "RSF", "unif"))

matplot(y = t(res1$surv_train[1:30,]), x = res1$time_points, type = "l")
print(res1$perf_train)
print(res1$max_time) # by default \code{max_time} is set to 2055 which is too large
# to have good predictions (training R2 with different weights are negative !)

# ------------------------------------------------
#   Training with estimation of test error
# ------------------------------------------------
set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res2 = cox_reg(y_var = "futime",
               delta_var = "delta",
               x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
               train = transplant_bis[train_lines,],
               test = transplant_bis[-train_lines,],
               types_w_ev = c("KM", "Cox", "RSF", "unif"))

print(res2$max_time) # default \code{max_time} has changed since train set
# is different

# train error is now positive but test error is still negative
print(res2$perf_train)
print(res2$perf_test)

# visualise the predictions
print(res2$pred_test[1:30])
matplot(y = t(res2$surv_test[1:30,]), x = res2$time_points, type = "l")

# ------------------------------------------------
#   Modify the \code{max_time} argument
# ------------------------------------------------
set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res3 = cox_reg(y_var = "futime",
               delta_var = "delta",
               x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
               train = transplant_bis[train_lines,],
               test = transplant_bis[-train_lines,],
               max_time = 600,
               types_w_ev = c("KM", "Cox", "RSF", "unif"))

print(res3$perf_train)
print(res3$perf_test) # test error is much better
print(res3$pred_test[1:30])
matplot(y = t(res3$surv_test[1:30,]), x = res3$time_points, type = "l")

# analyse the weights used for "weighted" criteria
print(res3$cens_rate) # rate of censoring taking into account \code{max_time}
print(head(res3$mat_w_test))
## ratio max(weights)/min(weights)
print(apply(X = res3$mat_w_test,
            MARGIN = 2,
            FUN = function(x){max(x[x != 0])/min(x[x != 0])}))
# ratios are low because the censoring rate is low

# in this case, it is not meaningful to to modify the
# \code{max_w_ev} argument since the maximum ratios
# between weights are around 2 and the test data has 197 rows.
# But in other situation it may be pertinent


# ------------------------------------------------
#   Use custom \code{phi} function
# ------------------------------------------------
g = function(x,a) abs(x-a)
set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res4 = cox_reg(y_var = "futime",
               delta_var = "delta",
               x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
               train = transplant_bis[train_lines,],
               test = transplant_bis[-train_lines,],
               phi = g,
               phi.args = list(a = 200), # set value for "a"
               max_time = 600,
               types_w_ev = c("KM", "Cox", "RSF", "unif"))

print(res4$perf_test)
print(res4$pred_test[1:30])



# ------------------------------------------------------------------------------------
#
#                             predict_cox_reg
#
# ------------------------------------------------------------------------------------

data("transplant", package = "survival")
transplant$delta = 1 * (transplant$event == "ltx") # create binary var
# which indicate censoring/non censoring

# keep only rows with no missing value
transplant_bis = transplant[stats::complete.cases(transplant),]


# ------------------------------------------------
#   Basic call to train a model
# ------------------------------------------------

set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res1 = cox_reg(y_var = "futime",
               delta_var = "delta",
               x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
               train = transplant_bis[train_lines,],
               types_w_ev = c("KM", "Cox", "RSF", "unif"))

# ------------------------------------------------
#   Predict on new data
# ------------------------------------------------

pred1 = predict_cox_reg(obj = res1,
                        newdata = transplant_bis[-train_lines,])
print(pred1$pred[1:30])





####################################################
####### autres tests à supprimer dans le futur

make_KM_weights(vect_y = veteran$time,
                vect_delta = veteran$status)

make_weights(data = veteran, y_name = "time", delta_name = "status", type = "Cox", max_ratio_weights = 1000,
             x_vars = setdiff(colnames(veteran),c("time","status")), cens_mod_obj = T)

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

predict_nodes2 = function (obj, newdata, na.action = stats::na.pass) {
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
    rpart:::pred.rpart(obj, as.matrix(newdata))
    predict(obj, newdata, type = "vector")
  }
  as.integer(row.names(obj$frame))[where]
}


library(rpart)
fit <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
predict_nodes(obj = fit, newdata = kyphosis, na.action = stats::na.pass)
predict_nodes2(obj = fit, newdata = kyphosis, na.action = stats::na.pass)

