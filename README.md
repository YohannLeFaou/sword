# sword : Survival Weighting fOr Regression of Duration

`sword` is an **experimental** package which implements the methodology that we study in [??]. The main function of sword is `sw_reg` and allows to use the Inverse Probability of Censoring Weighting (IPCW) method to fit random forest and gam model to &Phi;(T) , where T is a right censored variable and &Phi; is a given function.

Different types of weights are available (Kaplan-Meier, Cox, Random Survival Forest (RSF)) and custom weights may be provided. `cox_reg` and `rsf_reg` are functions which implement the 2 benchmarks based on Cox model and RSF model studied in our article. The package also provides different measures of performance regarding the quality of fit of the model.

`sword` can be used to get a quick insight about the precision of the IPCW method for the estimation of E[&Phi;(T)|X] when X is a vector of covariates, and to compare its performances with alternative models.

Technically, sword is a wrapper for statistical algorithms provided in the R packages [randomForestSRC](https://cran.r-project.org/web/packages/randomForestSRC/index.html), [rpart](https://cran.r-project.org/web/packages/rpart/index.html), [surv](https://cran.r-project.org/web/packages/surv/index.html) and [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html).

# Installation

sword is not available on CRAN and it can only be installed via devtools :

```
# First, install "devtools" :
install.packages("devtools")

# install "sword" :
devtools::install_github("YohannLeFaou/sword")
```

# Documentation

The [vignette of the sword package](https://cdn.rawgit.com/YohannLeFaou/sword/34d395ff/inst/doc/sword_vignette.html) explains the principle of the method and illustrates through an example how the package can be used.
