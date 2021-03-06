


# q_gamma() function -- inner loop of the variational BVS


# hyperparameters       : theta = (sigma, sigma * sa, pi), updated via MLE
# variational parameters: phi   = (alpha, mu, s), updated via cavi


## input:
#          X       : (n x p) design matrix, n is # of obs, p is # of variables
#          y       : (n x 1) response vector
#          sigma   : (1 x 1) variance of the RESIDUAL
#          sa      : (1 x 1) sigma * sa = prior var of REGRESSION COEFFICIENTS
#          logodds : (p x 1) prior log-odds of inclusion for each variable 
#          alpha   : (p x 1) var approx of P(gamma_k = 1 | X, y , theta)
#          mu      : (p x 1) mean of coeff_i given that it's in model, mu0[i]

## output:
#          logw    : (iter x 1) variational estimate of the marginal log-like
#          err     : (iter x 1) diff between alphas between iterations
#          sigma   : (1 x 1)    updated (if asked) variance of residual
#          sa      : (1 x 1)    updated (if asked) variance of reg. coeff
#          alpha   : (p x 1)    updated alpha (see above input description)
#          mu      : (p x 1)    updated mu    (see above input description)
#          s       : (s x 1)    updated posterior variance


## notes on output:
## q_gamma is called once per iteration of the outer loop, outputting an alpha
## per iteration. These alpha's (each is a p-dimensional vector) are weighted
## to get a Monte Carlo estimate of the posterior inclusion probability of each
## of the regression coefficients


## quantities that are calculated within the code:
#          iter    : number of iterations to achieve convergence
#          Xr      : efficient way to store X * E(beta | gamma = 1)


# to do: figure out how to modify/interpret roles of n0, sa0

q_gamma = function(X, y, sigma, sa, logodds, alpha, mu,
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
		logw = varLB(Xr, d, y, sigma, alpha, mu, s, logodds, sa) # calculate.R


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
		## these still need to be derived -- why are these updated after VLB
		if (update.sigma) {
			sigma = l2(y - Xr)^2 + sum(d * var_ss(alpha, mu, s)) + 
			        sum(alpha * (s + mu^2) / sa) / (n + sum(alpha))
			s     = sa * sigma / (d * sa + 1)
		}

		if (update.sa) {
			sa_numer = sa0 * n0 + sum(alpha * (s + mu^2))
			sa       = sa_numer / (n0 + sigma * sum(alpha))
			s        = sa * sigma / (d * sa + 1)
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
	out = list(vlb = logw[iter], logw = logw[1:iter], err = err[1:iter], 
		     sigma = sigma, sa = sa, alpha = alpha, mu = mu, s = s, iter = iter)

	return(out)

} # end q_gamma() function










