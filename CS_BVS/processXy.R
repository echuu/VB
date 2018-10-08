# Adjust variables X and continuous outcome Y so that the linear
# effects of the covariates Z are removed. This is equivalent to
# integrating out the regression coefficients corresponding to the
# covariates with respect to an improper, uniform prior; see Chipman,
# George and McCulloch, "The Practical Implementation of Bayesian
# Model Selection," 2001. It is assumed that the first column of Z is
# the intercept; that is, a column of ones.

processXy <- function (X, y) {

  
  X = scale(X,center = TRUE,scale = FALSE)
  y = y - mean(y)

  return(list(X = X,y = y))
} # end of processXy() function

# end of processXy.R file
