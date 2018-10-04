# ----------------------------------------------------------------------
# This function implements one iteration of the "outer loop".
outerloop <- function (X, Z, y, family, weights, resid.vcov, SZy, SZX, sigma,
                       sa, logodds, alpha, mu, eta, update.order, tol, maxiter,
                       verbose, outer.iter, update.sigma, update.sa,
                       optimize.eta, n0, sa0) {
  p <- ncol(X)
  if (length(logodds) == 1)
    logodds <- rep(logodds,p)

  # Note that we need to multiply the prior log-odds by log(10),
  # because varbvsnorm, varbvsbin and varbvsbinz define the prior
  # log-odds using the natural logarithm (base e).
  if (family == "gaussian") {

      # Optimize the variational lower bound for the Bayesian variable
      # selection model.
      out <- varbvsnorm(X,y,sigma,sa,log(10)*logodds,alpha,mu,update.order,
                        tol,maxiter,verbose,outer.iter,update.sigma,update.sa,
                        n0,sa0)

      # Adjust the variational lower bound to account for integral over
      # the regression coefficients corresponding to the covariates.
      out$logw <- out$logw - determinant(crossprod(Z),logarithm = TRUE)$modulus/2
      
      # Compute the posterior mean estimate of the regression
      # coefficients for the covariates under the current variational
      # approximation.
      # this calculation is slightly modified version of (9) in C&S (2012)
      out$mu.cov <- c(with(out,SZy - SZX %*% (alpha*mu)))
  } 

  numiter  <- length(out$logw)
  out$logw <- out$logw[numiter]
  return(out)
}
