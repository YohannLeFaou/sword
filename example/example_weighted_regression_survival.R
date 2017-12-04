
# ------------------------------------------------------------------------------------
#
#                             weighted_regression_survival
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
res1 = weighted_regression_survival(y_var = "futime",
                                    delta_var = "delta",
                                    x_vars = setdiff(colnames(transplant_bis),
                                                     c("futime", "delta", "event")),
                                    data_train = transplant_bis
)
# parameters set by default
res1$type_weights
res1$type_regression
res1$max_ratio_weights_model
res1$max_ratio_weights_eval
res1$mode_w_RF # 1 corresponds to wRF1 in [Gerb. et al.]


# train errors
res1$list_criteria_train

# ------------------------------------------------
#   Training with estimation of test error
# ------------------------------------------------
set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res2 = weighted_regression_survival(y_var = "futime",
                                    delta_var = "delta",
                                    x_vars = setdiff(colnames(transplant_bis),
                                                     c("futime", "delta", "event")),
                                    data_train = transplant_bis[train_lines,],
                                    data_test = transplant_bis[-train_lines,],
                                    types_weights_eval = c("KM", "Cox", "RSF", "unif"))

print(res2$max_time) # default \code{max_time} has changed since train set
# is different

# there is a uge overfitting in terms of quadratic errors
print(res2$list_criteria_train)
print(res2$list_criteria_test)

# default parameters for the random forest are
res2$weighted_RF_object$ntree
res2$weighted_RF_object$mtry
res2$weighted_RF_object$nodesize
res2$weighted_RF_object$nodedepth # means there is no depth limit


# -----------------------------------------------------
#   Modify the \code{max_time} argument & look for
#       the best model under this setting
# -----------------------------------------------------

set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res30 = weighted_regression_survival(y_var = "futime",
                                    delta_var = "delta",
                                    x_vars = setdiff(colnames(transplant_bis),
                                                     c("futime", "delta", "event")),
                                    data_train = transplant_bis[train_lines,],
                                    data_test = transplant_bis[-train_lines,],
                                    type_weights = "KM", # default value
                                    max_time = 600, # we set \code{max_time} to 600
                                    types_weights_eval = c("KM", "Cox", "RSF", "unif"))
print(res30$list_criteria_test)

# are the other types of weights giving better results ?
## Cox weights
res31 = weighted_regression_survival(y_var = "futime",
                                    delta_var = "delta",
                                    x_vars = setdiff(colnames(transplant_bis),
                                                     c("futime", "delta", "event")),
                                    data_train = transplant_bis[train_lines,],
                                    data_test = transplant_bis[-train_lines,],
                                    type_weights = "Cox",
                                    max_time = 600, # we set \code{max_time} to 600
                                    types_weights_eval = c("KM", "Cox", "RSF", "unif"))
print(res31$list_criteria_test) # slight improvment compared with weights KM

## RSF weights
res32 = weighted_regression_survival(y_var = "futime",
                                     delta_var = "delta",
                                     x_vars = setdiff(colnames(transplant_bis),
                                                      c("futime", "delta", "event")),
                                     data_train = transplant_bis[train_lines,],
                                     data_test = transplant_bis[-train_lines,],
                                     type_weights = "RSF",
                                     max_time = 600, # we set \code{max_time} to 600
                                     types_weights_eval = c("KM", "Cox", "RSF", "unif"))
print(res32$list_criteria_test)

## unif weights
res33 = weighted_regression_survival(y_var = "futime",
                                     delta_var = "delta",
                                     x_vars = setdiff(colnames(transplant_bis),
                                                      c("futime", "delta", "event")),
                                     data_train = transplant_bis[train_lines,],
                                     data_test = transplant_bis[-train_lines,],
                                     type_weights = "unif",
                                     max_time = 600, # we set \code{max_time} to 600
                                     types_weights_eval = c("KM", "Cox", "RSF", "unif"))
print(res33$list_criteria_test)

# In terms of quadratic, the best weights are the Cox weights
# remark : in this example there is not a big difference between the "unif" weights
# and the other weights because there is little censoring in the data :
res33$censoring_rate_with_threshold

# -----------------------------------------------------------------
#     Try wRF2  and wRF3 (both are obtained with
#               \code{mode_w_RF = 2})
# -----------------------------------------------------------------

res40 = weighted_regression_survival(y_var = "futime",
                                     delta_var = "delta",
                                     x_vars = setdiff(colnames(transplant_bis),
                                                      c("futime", "delta", "event")),
                                     data_train = transplant_bis[train_lines,],
                                     data_test = transplant_bis[-train_lines,],
                                     type_weights = "Cox",
                                     max_time = 600, # we set \code{max_time} to 600
                                     types_weights_eval = c("KM", "Cox", "RSF", "unif"),
                                     mode_w_RF = 2)
print(res40$list_criteria_test) # wRF2 : not as good as wRF1
print(res40$list_criteria_test_KMloc) # wRF3 : worse than wRF2


# -------------------------------------------------------
#               Try a GAM model
# -------------------------------------------------------

## GLM with Cox weights
res5 = weighted_regression_survival(y_var = "futime",
                                     delta_var = "delta",
                                     x_vars = setdiff(colnames(transplant_bis),
                                                      c("futime", "delta", "event")),
                                     data_train = transplant_bis[train_lines,],
                                     data_test = transplant_bis[-train_lines,],
                                     type_weights = "Cox",
                                     max_time = 600, # we set \code{max_time} to 600
                                     types_weights_eval = c("KM", "Cox", "RSF", "unif"),
                                    type_regression = "gam")
print(res5$list_criteria_test) # not as good as random forest


# ------------------------------------------------------------
#     Analyse the weights used for "weighted" criteria
# ------------------------------------------------------------

print(res31$censoring_rate_with_threshold) # rate of censoring taking into account \code{max_time}
print(head(res31$mat_weights_test))
## ratio max(weights)/min(weights)
print(apply(X = res31$mat_weights_test,
            MARGIN = 2,
            FUN = function(x){max(x[x != 0])/min(x[x != 0])}))
# ratios are low because the censoring rate is low

# in this case, it is not meaningful to to modify the
# \code{max_ratio_weights_eval} argument since the maximum ratios
# between weights are around 2 and the test data has 197 rows.
# But in other situation it may be pertinent


# ------------------------------------------------
#   Use custom \code{phi} function
# ------------------------------------------------
g = function(x,a) abs(x-a)
set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
res6 = weighted_regression_survival(y_var = "futime",
                                    delta_var = "delta",
                                    x_vars = setdiff(colnames(transplant_bis),
                                                     c("futime", "delta", "event")),
                                    data_train = transplant_bis[train_lines,],
                                    data_test = transplant_bis[-train_lines,],
                                    phi = g,
                                    phi.args = list(a = 200),
                                    type_weights = "Cox",
                                    max_time = 600, # we set \code{max_time} to 600
                                    types_weights_eval = c("KM", "Cox", "RSF", "unif"))
print(res6$list_criteria_test) # slight improvment compared with weights KM



# ------------------------------------------------------------------------------------
#
#                             predict_weighted_regression_survival
#
# ------------------------------------------------------------------------------------

# ------------------------------------------------
#   Load "transplant" data
# ------------------------------------------------
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
res = weighted_regression_survival(y_var = "futime",
                                    delta_var = "delta",
                                    x_vars = setdiff(colnames(transplant_bis),
                                                     c("futime", "delta", "event")),
                                    data_train = transplant_bis[train_lines,],
                                    types_weights_eval = c("KM", "Cox", "RSF", "unif"),
                                    mode_w_RF = 2)
# ------------------------------------------------
#   Predict on new data
# ------------------------------------------------

pred = predict_weighted_regression_survival(object = res,
                                            newdata = transplant_bis[-train_lines,])

