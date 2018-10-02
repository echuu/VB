
# bvs.R 

## Bayesian Variable Selection via the method described in Carbonetto & Stephens
## Each iteration gives an alpha, which are used to generate an MCMC estimate of
## the posterior inclusion probability of the regression coefficients


bvs = function(X, y, sigma, sa, logodds,
	           update.sigma, update.sa, B = 100, sa0 = 1, n0 = 10,
	           tol = 1e-4, maxiter = 1e4) {

} # end of bvs() functio







