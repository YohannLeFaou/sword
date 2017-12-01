
# ------------------------------------------------------------------------------------
#
#                             NormalizedGini
#
# ------------------------------------------------------------------------------------
set.seed(17)
x = runif(1000)
y = rnorm(mean = x, sd = 1, n = 1000)

NormalizedGini(solutions = x,
               predictions = y)

NormalizedGini(solutions = x,
               predictions = y,
               weights = x)

# be carefull
NormalizedGini(solutions = x,
               predictions = y,
               weights = abs(y))
# when the weights depends on \code{predictions}, \code{NormalizedGini}
# may gives strange results

# ------------------------------------------------------------------------------------
#
#                             eval_weighted_criteria
#
# ------------------------------------------------------------------------------------
set.seed(17)
x = runif(1000)
y = rnorm(mean = x, sd = 0.2, n = 1000)

eval_weighted_criteria(predictions = y,
                       solutions = x)

