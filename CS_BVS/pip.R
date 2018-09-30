



# mu[i], s[i] are the mean/variance of the coefficient given that it's in model


pip = function(X, y, sigma, sa, logodds, alpha, mu, update.order,
	           tol = 1e-4, maxiter = 1e4, verbose = TRUE, 
	           outer.iter = NULL, update.sigma = TRUE,
	           update.sa = TRUE, n0 = 10, sa0 = 1) {


	n = nrow(X)
	p = ncol(X)


	# pre-computations:
	xy = c(y %*% X)
	diag = diag(t(X) %*% X)
	Xr = c(X %*% (alpha * mu))

	s = sa * sigma / (sa * d + 1) # initial value for Var(beta_k | gamma_k = 1)

	logw = rep(0, maxiter) # variational estimate of the marginal log-like
	err  = rep(0, maxiter) # max diff b/w approimxate posterior prob (alpha)
	                       # at successive iterations, also used to check conv.


	# main loop
	for (i in 1:maxiter) {
		
		# variational parameters
		alpha0 = alpha
		mu0    = mu      # posterior mean of coefficients
		s0     = s       # posterior variance of coefficients

		# hyperparameters
		sigma0 = sigma   # variance of residuals
		sa.old = sa      # prior variance of coefficients

		# compute variational lower bound
		logw = varLB(Xr, d, y, sigma, alpha, mu, s, logodds, sa)

		# update variational parameters via CAVI, phi = (alpha, mu, s)



		# carbonetto code doe



	} # end for() loop








} # end pip() function










