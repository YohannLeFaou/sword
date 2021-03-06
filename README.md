# sword : Survival Weighting fOr Regression of Duration

`sword` is a R package which implements the methodology that we study in our [article](https://hal.archives-ouvertes.fr/hal-01807623/file/Article___The_impact_of_churn_on_prospect_value_in_health_insurance__evaluation_using_a_random_forest_under_random_censoring%20%281%29.pdf). The main function of sword is `sw_reg`. It allows to use the Inverse Probability of Censoring Weighting (IPCW) method to fit random forest and gam model to &Phi;(T) , where T is a right censored variable and &Phi; is a given function.

Different types of weights are available (Kaplan-Meier, Cox, Random Survival Forest (RSF)) and custom weights may be provided. `cox_reg`,  `rsf_reg`, `rrt_reg` and `rlt_reg` are functions which implement benchmarks models based on respectively: the Cox model, the Random Survival Forest algorithm, the Relative risk tree model, and the Reinforcement Learning Tree algorithm, all studied in our article. The package also provides different performance measures of the quality of fit of the models.

`sword` can be used to get a quick insight about the precision of the IPCW method for the estimation of E[&Phi;(T)|X] when X is a vector of covariates, and to compare its performances with alternative models.

Technically, sword is a wrapper for statistical algorithms provided in the R packages [randomForestSRC](https://cran.r-project.org/web/packages/randomForestSRC/index.html), [rpart](https://cran.r-project.org/web/packages/rpart/index.html), [surv](https://cran.r-project.org/web/packages/surv/index.html), [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html), and [RLT](https://cran.r-project.org/web/packages/RLT/index.html).

# Installation

sword is not available on CRAN and it can only be installed via devtools :

```
# First, install "devtools" :
install.packages("devtools")

# install "sword" :
devtools::install_github("YohannLeFaou/sword")
```

# Documentation

The [vignette of the sword package](https://rawgit.com/YohannLeFaou/sword/master/inst/doc/sword_vignette.html) explains the principle of the method and illustrates through an example how the package can be used.
