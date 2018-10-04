# Adjust variables X and continuous outcome Y so that the linear
# effects of the covariates Z are removed. This is equivalent to
# integrating out the regression coefficients corresponding to the
# covariates with respect to an improper, uniform prior; see Chipman,
# George and McCulloch, "The Practical Implementation of Bayesian
# Model Selection," 2001. It is assumed that the first column of Z is
# the intercept; that is, a column of ones.

processXy <- function (X, y, Z) {

  # Here I compute two quantities that are used here to remove linear
  # effects of the covariates (Z) on X and y, and later on to
  # efficiently compute estimates of the regression coefficients for
  # the covariates.
  A   = t(Z) %*% Z
  SZy = as.vector(solve(A,c(y %*% Z)))
  SZX = as.matrix(solve(A,t(Z) %*% X))
  
  X = scale(X,center = TRUE,scale = FALSE)
  y = y - mean(y)

  #if (ncol(Z) == 1) {
  #  X <- scale(X,center = TRUE,scale = FALSE)
  #  y <- y - mean(y)
  #} else {

    # The equivalent expressions in MATLAB are  
    #
    #   y = y - Z*((Z'*Z)\(Z'*y))
    #   X = X - Z*((Z'*Z)\(Z'*X))  
    #
    # This should give the same result as centering the columns of X
    # and subtracting the mean from y when we have only one
    # covariate, the intercept.
    #y <- y - c(Z %*% SZy)
    #X <- X - Z %*% SZX
  #}

  return(list(X = X,y = y,SZy = SZy,SZX = SZX))
} # end of processXy() function

# end of processXy.R file
