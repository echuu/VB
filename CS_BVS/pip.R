



# mu[i], s[i] are the mean/variance of the coefficient given that it's in model


## carbonetto documentation seems to have different formulation of sa
## from the prior variance described in the paper -- check results

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
		
		#### ----   store variational parameters from prev iteration   ---- ####
		# ----------------------------------------------------------------------
		alpha0 = alpha
		mu0    = mu      # posterior mean of coefficients
		s0     = s       # posterior variance of coefficients


		#### ----       store hyperparameters from prev iteration      ---- ####
		# ----------------------------------------------------------------------
		sigma0 = sigma   # variance of residuals
		sa.old = sa      # prior variance of coefficients


		#### ----           compute variational lower bound            ---- ####
		# ----------------------------------------------------------------------
		logw = varLB(Xr, d, y, sigma, alpha, mu, s, logodds, sa)


		#### ----         perform one cycle of cavi updates            ---- ####
		# variational parameters, phi = (alpha, mu, s), are updated via
		# call to vbUpdate(), after one cycle of updates, we update values for
		# alpha, mu, Xr so that we can recompute the variational lower bound
		# ----------------------------------------------------------------------
		cavi_update = vbUpdate(X, sigma, sa, logodds, xy, d, alpha, mu, Xr)

		alpha = cavi_update$alpha    # update alpha
		mu    = cavi_update$mu       # update mu
		Xr    = cavi_update$Xr       # update Xr


		#### ----        recompute the variational lower bound         ---- ####
		# ----------------------------------------------------------------------		
		logw[i] = varLB(Xr, d, y, sigma, alpha, mu, s, logodds, sa)


		#### ----                update hyperparameters                ---- ####
		# ----------------------------------------------------------------------

		if (update.sigma) {
			
		}

		if (update.sa) {
			
		}


		#### ----                  check convergence                   ---- ####
		# ----------------------------------------------------------------------



	} # end for() loop


	#### ----                    prepare output                        ---- ####
	# --------------------------------------------------------------------------



} # end pip() function










