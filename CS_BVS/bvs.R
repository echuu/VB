
# bvs.R 

## Bayesian Variable Selection via the method described in Carbonetto & Stephens
## Each iteration gives an alpha, which are used to generate an MCMC estimate of
## the posterior inclusion probability of the regression coefficients


## input:
#          X       : (n x p) design matrix, n is # of obs, p is # of variables
#          sigma   : (1 x 1) variance of the RESIDUAL
#          sa      : (1 x 1) PRIOR variance of regression coefficients
#          logodds : (p x 1) prior log-odds of inclusion for each variable 

library(dplyr)


bvs = function(X, y, sigma, sa, logodds,
	           update.sigma, update.sa, B = 100, sa0 = 1, n0 = 10,
	           tol = 1e-4, maxiter = 1e4) {

    n = nrow(X)
	p = ncol(X)

	# carbonetto code allows for 
		# (1) variance of residual
		# (2) prior variance of regression coefficients 
		# (3) logodds 
    # to be missing in looks for candidate settings, but for now, we require
    # these parameters to be present

	#### ----                      check input                         ---- ####
	# --------------------------------------------------------------------------
	checkInitialInput(X, y, sigma, sa, logodds, n, p)

	y = c(as.double(y))   # y must be in this form and not stored as matrix,
	                      # else processXy() function will complain about dim.

	#### ----                  preprocess variables                    ---- ####
	# --------------------------------------------------------------------------

	# create intercept
	Z = matrix(1, n, 1)

    # Adjust the inputs X and y so that the linear effects of the
    # covariates (Z) are removed. This is equivalent to integrating
    # out the regression coefficients corresponding to the covariates
    # with respect to an improper, uniform prior.	

    new_Xy = processXy(X, y, Z)

    # for initial applications, Z will be the intercept, so the resulting X will
    # just be the centered version, and y = y - mean(y)
    X   = new_Xy$X
    y   = new_Xy$y 
    
    SZX = new_Xy$SZX    # only used if there are additional covariates
	SZy = new_Xy$SZy    # only used if there are additional covariates

    rm(new_Xy)

	#### ----            initialize variational parameters             ---- ####
	# --------------------------------------------------------------------------
	alpha = runif_mat(p, B)    # (p x B) estimates of each iter stored col-wise
	mu    = runif_mat(p, B) 



	# carbonetto does something here with a Z matrix, integrating out coeffs
	
	#### ----                  preprocess variables                    ---- ####
	# --------------------------------------------------------------------------	# initialize output variables
	logw = rep(0, B)                   # var. estimate of marginal log-like
	s    = matrix(0, p, B)             # variance of reg coeffs, stored col-wise
	# mu.cov = matrix(0, ncol(Z), B)   # posterior mean estimates
	    # how is mu.cov different from mu matrix declared in previous section




    #### ----                       outer loop                         ---- ####
	# --------------------------------------------------------------------------
	for (i in 1:B) {

		# inner loop -- optimize var lower bound
		cavi = q_gamma(X, y, sigma, sa, logodds, alpha, mu) 


		# compute posterior mean estimate of the regression coefficients
		# for the covariates under the current variational approximation
		# add posterior mean as a new parameter to output
		post_mean = c()
		cavi = out %>% mutate(post_mean = post_mean)


		# store updates for i-th iteration
		logw[i]    = cavi$vlb
		sigma[i]   = cavi$sigma
		sa[i]      = cavi$sa
		mu.cov[,i] = cavi$post_mean # posterior mean ?
		alpha[,i]  = cavi$alpha
		mu[,i]     = cavi$mu        # in q_gamma, this is computed as post. mean
		s[,i]      = cavi$s

	} # end for()



	# prepare output


} # end of bvs() function







