## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(sword)
library(ggplot2)
library(reshape2)

## ------------------------------------------------------------------------
data("transplant")
head(transplant)

## ---- fig.show='hold'----------------------------------------------------
transplant$delta = 1 * (transplant$event == "ltx") # create binary variable which indicate censoring/non censoring

## ---- fig.show='hold'----------------------------------------------------
# compute the number of missing values for each column
apply(transplant, MARGIN = 2, FUN = function(x){sum(is.na(x))})

# keep only rows with no missing value
transplant_bis = transplant[complete.cases(transplant),]

## ---- fig.width=7, fig.height=4------------------------------------------
# plot the survival curve of the waiting time until transplant
KM = survfit(formula = Surv(time = futime, event = delta) ~ 1,
             data = transplant_bis)
plot(KM, 
     main = "Survival Curve of the waiting time before liver transplant", 
     ylab = "survival prob.", xlab = "time (days)")

## ------------------------------------------------------------------------
nrow(transplant_bis) # number of observations

set.seed(17)
train_lines = sample(1:nrow(transplant_bis), 600)
train = transplant_bis[train_lines,] # train set
test = transplant_bis[-train_lines,] # test set

## ------------------------------------------------------------------------
res1 = sw_reg(y_var = "futime", 
              delta_var = "delta",
              x_vars = c("age", "sex", "abo", "year"),
              train = train,
              test = test,
              type_w = "KM",
              phi = function(x){(x > 365) * 1},
              max_time = 366,
              types_w_ev = c("KM", "Cox", "RSF", "unif"))

## ------------------------------------------------------------------------
print(sum(train$delta == 0) / nrow(train)) # rate of censoring on the initial data
print(head(res1$cens_rate)) # rate of censoring after thresholding with `max_time`

## ------------------------------------------------------------------------
print(res1$w_mod_train[1:30])

## ------------------------------------------------------------------------
print(res1$max_w_mod) # value passed to the sw_reg function (default value = 20)
print(max(res1$w_mod_train) / min(res1$w_mod_train[res1$w_mod_train > 0])) # maximum ratio among the train weights
print(res1$n_w_mod_modif_train) # number of weights modified due to `max_w_mod`

## ------------------------------------------------------------------------
print(head(res1$pred_train))
print(head(res1$pred_test))

## ------------------------------------------------------------------------
pred_test = predict_sw_reg(obj = res1, newdata = test)
print(pred_test$pred[1:30])

## ------------------------------------------------------------------------
print(res1$perf_test$criteria_weighted) # test mse and R2

## ------------------------------------------------------------------------
print(res1$perf_test$concordance)

## ------------------------------------------------------------------------
print(res1$perf_train)

## ------------------------------------------------------------------------
# weights used for evaluation
print(head(res1$mat_w_train))
print(head(res1$mat_w_test))

## ------------------------------------------------------------------------
print(res1$max_w_ev) # value passed to the sw_reg function (default value = 1000)
print(res1$n_w_ev_modif_test) # number of test weights modified because of the threshold
print(res1$n_w_ev_modif_train) # number of train weights modified because of the threshold

## ------------------------------------------------------------------------
print(res1$sum_w_train) # sum of the train weights before reprocessing
print(res1$sum_w_test) # sum fo the test weighs before reprocessing

## ------------------------------------------------------------------------
res2 = sw_reg(y_var = "futime",
              delta_var = "delta",
              x_vars = c("age", "sex", "abo", "year"),
              train = train,
              test = test,
              type_w = "KM",
              phi = function(x){(x > 365) * 1},
              max_time = 366,
              types_w_ev = c("KM", "Cox", "RSF", "unif"),
              mode_sw_RF = 2)

## ------------------------------------------------------------------------
print(res2$perf_test)
print(res2$perf_test_KMloc)

## ------------------------------------------------------------------------
res2 = sw_reg(y_var = "futime", 
              delta_var = "delta",
              x_vars = c("age", "sex", "abo", "year"),
              train = train,
              test = test,
              type_reg = "gam",
              type_w = "Cox",
              phi = function(x){(x > 365) * 1},
              max_time = 366,
              types_w_ev = c("KM", "Cox", "RSF", "unif"),
              family = binomial(link = "logit"))

## ------------------------------------------------------------------------
summary(res2$sw_gam_obj)

## ------------------------------------------------------------------------
print(res2$pred_test[1:20])
print(res2$perf_test)

## ------------------------------------------------------------------------
res11 = sw_reg(y_var = "futime", 
               delta_var = "delta",
               x_vars = c("age", "sex", "abo", "year"),
               train = train,
               test = test,
               type_w = "KM",
               phi = function(x){(x > 365) * 1},
               max_time = 366,
               types_w_ev = c("KM", "Cox", "RSF", "unif"),
               proximity = T)

## ------------------------------------------------------------------------
print(res11$sw_RF_obj$proximity[1:5,1:5]) # matrix of the proximities between the first 5 obs. of the train set
print(dim(res11$sw_RF_obj$proximity)) # dimension of the proximity matrix

## ------------------------------------------------------------------------
res2 = rsf_reg(y_var = "futime",
               delta_var = "delta",
               x_vars = c("age", "sex", "abo", "year"),
               train = transplant_bis[train_lines,],
               test = transplant_bis[-train_lines,],
               phi = function(x){(x > 365) * 1},
               max_time = 366,
               types_w_ev = c("KM", "Cox", "RSF", "unif"))

## ---- fig.width=7--------------------------------------------------------
data_surv_test = cbind(melt(t(res2$surv_test[1:10,])), time = res2$time_points)
ggplot(data = data_surv_test, 
       aes(x = time, y = value, group = factor(Var2), color = factor(Var2))) +
  geom_line() +
  theme(legend.position = "bottom") +
  ggtitle("Prédiction des courbes de survie des 10 premiers individus test (RSF)")

## ------------------------------------------------------------------------
res2$pred_test[1:10]

## ------------------------------------------------------------------------
print(res2$perf_test)

## ------------------------------------------------------------------------
res3 = cox_reg(y_var = "futime",
               delta_var = "delta",
               x_vars = c("age", "sex", "abo", "year"),
               train = transplant_bis[train_lines,],
               test = transplant_bis[-train_lines,],
               phi = function(x){(x > 365) * 1},
               max_time = 366,
               types_w_ev = c("KM", "Cox", "RSF", "unif")
               )

## ----  fig.width=7-------------------------------------------------------
data_surv_test2 = cbind(melt(t(res3$surv_test[1:10,])), time = res3$time_points)
ggplot(data = data_surv_test2, 
       aes(x = time, y = value, group = factor(Var2), color = factor(Var2))) +
  geom_line() +
  theme(legend.position = "bottom") +
  ggtitle("Prédiction des courbes de survie des 10 premiers individus test (Cox)")

## ------------------------------------------------------------------------
print(res3$perf_test)

## ------------------------------------------------------------------------
# là j'ai calculé le R2 du modele qui est la moyenne entre weighted rf et RSF, mai visiblement pas ouf (trop de variance dans les résultats de tte façon)
mean = sum( res1$mat_w_test[,"KM"] * res1$test$phi)
R2 = 1 - sum( res1$mat_w_test[,"KM"] * ((res1$pred_test + res2$pred_test)/2 - res1$test$phi )^2) / 
  sum( res1$mat_w_test[,"KM"] * (res1$test$phi - mean )^2)

print(mean)
print(R2)

