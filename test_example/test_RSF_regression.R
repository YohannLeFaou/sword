
# ------------------------------------------------------------------------------------
#
#                             RSF_regression
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

res1 = RSF_regression(y_var = "futime",
                      delta_var = "delta",
                      x_vars = setdiff(colnames(transplant_bis),
                                       c("futime", "delta", "event")),
                      data_train = transplant_bis,
                      types_weights_eval = c("KM", "Cox", "RSF", "unif"))

# by default, (main) parameters used for the random survival forest are :
print(res1$RSF_object$mtry)
print(res1$RSF_object$nodesize)
print(res1$RSF_object$nodedepth) # "-1" = no depth limitation

# visualise the train predictions
matplot(y = t(res1$survival_train[1:30,]), x = res1$time_points, type = "l")
print(res1$list_criteria_train)

print(res1$max_time) # by default \code{max_time} is set to 2055 which is very large
# given the outlook of the survival function of \eqn{T}. Train errors may be
# overfitted

# ------------------------------------------------
#   Training with estimation of test error
# ------------------------------------------------
set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res2 = RSF_regression(y_var = "futime",
                      delta_var = "delta",
                      x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
                      data_train = transplant_bis[train_lines,],
                      data_test = transplant_bis[-train_lines,],
                      types_weights_eval = c("KM", "Cox", "RSF", "unif"))

print(res2$max_time) # default \code{max_time} has changed since the train set
# is different

# train error is positive but test error is negative
print(res2$list_criteria_train)
print(res2$list_criteria_test) # weighted criterias show there is a lot of overfitting

# ------------------------------------------------
#   Modify the \code{max_time} argument
# ------------------------------------------------
set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res3 = RSF_regression(y_var = "futime",
                      delta_var = "delta",
                      x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
                      data_train = transplant_bis[train_lines,],
                      data_test = transplant_bis[-train_lines,],
                      max_time = 600,
                      types_weights_eval = c("KM", "Cox", "RSF", "unif"))

print(res3$list_criteria_train)
print(res3$list_criteria_test) # test error is much better

# visualise the predictions
print(res3$predicted_test[1:30])
matplot(y = t(res3$survival_test[1:30,]), x = res3$time_points, type = "l")

# analyse the weights used for "weighted" criteria
print(res3$censoring_rate_with_threshold) # rate of censoring taking into account \code{max_time}
print(head(res3$mat_weights_test))
## ratio max(weights)/min(weights)
print(apply(X = res3$mat_weights_test,
            MARGIN = 2,
            FUN = function(x){max(x[x != 0])/min(x[x != 0])}))
# ratios are low because the censoring rate is low

# in this case, it is not meaningful to to modify the
# \code{max_ratio_weights_eval} argument since the maximum ratios
# between weights are around 2 and the test data has 197 rows.
# But in other situation it may be pertinent

# ----------------------------------------------------------------
#   Change the parameter for \code{\link[randomForestSRC]{rfsrc}}
# ----------------------------------------------------------------
set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res4 = RSF_regression(y_var = "futime",
                      delta_var = "delta",
                      x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
                      data_train = transplant_bis[train_lines,],
                      data_test = transplant_bis[-train_lines,],
                      max_time = 600,
                      types_weights_eval = c("KM", "Cox", "RSF"),
                      minleaf = 5 # change \code{nodesize} for the inner call to \code{\link[randomForestSRC]{rfsrc}}
                      )

print(res4$list_criteria_test) # slight amelioration compared


# ------------------------------------------------
#   Use custom \code{phi} function
# ------------------------------------------------
g = function(x,a) abs(x-a)
set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res5 = RSF_regression(y_var = "futime",
                      delta_var = "delta",
                      x_vars = setdiff(colnames(transplant_bis),c("futime", "delta", "event")),
                      data_train = transplant_bis[train_lines,],
                      data_test = transplant_bis[-train_lines,],
                      phi = g,
                      phi.args = list(a = 200), # set value for "a"
                      max_time = 600,
                      types_weights_eval = c("KM", "Cox", "RSF", "unif"))

print(res5$list_criteria_test)
print(res5$predicted_test[1:30])

