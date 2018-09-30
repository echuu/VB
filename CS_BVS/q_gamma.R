


# q_gamma() function -- inner loop of the variational BVS


# mu[i], s[i] are the mean/variance of the coefficient given that it's in model


## carbonetto documentation seems to have different formulation of sa
## from the prior variance described in the paper -- check results


## input:
#          X       : (n x p) design matrix, n is # of obs, p is # of variables
#          y       : (n x 1) response vector
#          sigma   : (1 x 1) variance of the RESIDUAL
#          sa      : (1 x 1) prior variance of REGRESSION COEFFICIENTS
#          logodds : (p x 1) prior log-odds of inclusion for each variable 
#          alpha   : (p x 1) curr param of var approx of PIP(i) = alpha0[i]
#          mu      : (p x 1) mean of coeff_i given that it's in model, mu0[i]



q_gamma = function(X, y, sigma, sa, logodds, alpha, mu, update.order,
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
	iter = 0               # number of iterations for convergence

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

		# check 1: if variational lower bound of prev iteration > updated 
		#          variational lower bound, use prev. estimates and EXIT
		VLB_CONVERGE = logw[i] < logw0

		## check 2: if the change in the error vector < tol, then EXIT 
		## update the change in alpha
		err[i]         = max(abs(alpha - alpha0))
		ERROR_CONVERGE = err[i] < tol

		if (VLB_CONVERGE) {

			logw[i]   = logw0
			err[i]    = 0
			sigma     = sigma0
			sa        = sa.old
			alpha     = alpha0
			mu        = mu0
			s         = s0
			iter      = i

			break
		}

		if (ERROR_CONVERGE) {
			iter = i
			break
		}

	} # end of main for() loop


	#### ----                    prepare output                        ---- ####
	# --------------------------------------------------------------------------
	out = list(logw = logw[1:iter], err = err[1:iter], sigma = sigma,
		       sa = sa, alpha = alpha, mu = mu, s = s)

	return(out)

} # end q_gamma() function










