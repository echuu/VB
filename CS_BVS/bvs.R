
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


	# scale the dimensions of the (hyper/variational)-parameters so that we can
	# index into them as vectors

	
	#### ----                  preprocess variables                    ---- ####
	# --------------------------------------------------------------------------	

	# initialize output variables
	logw = rep(0, B)                   # var. estimate of marginal log-like
	s    = matrix(0, p, B)             # variance of reg coeffs, stored col-wise


    #### ----                       outer loop                         ---- ####
	# --------------------------------------------------------------------------
	
	# (1) If initialization of the parameters is not provided in the input, 
	#     we iterate through the variational updates to find a better 
	#     initialization for the variational parameters. The values we use to
	#     run CAVI for this chunk are generated randomly when we defined above
	for (i in 1:B) {



		## find the iteration with parameters that maximized the VLB
		## store these parameters in alpha, mu, sigma (*), sa (*)
		## (*): if update.sigma == TRUE

	} # end of initilization for()



	# (2) Using the initialization chosen in step (1), we run CAVI.
	#     
	for (i in 1:B) {

		# inner loop -- optimize var lower bound
		if (length(logodds) == 1) {
			logodds = rep(logodds, p)
		}

		## details of these updates still need to be figured out..
		## seems that if we have initialize.params == TRUE, then all starting
		## values of alpha[,i], mu[,i] will all be the same, which would result
		## in the same updates stored in each iteration, which would defeat the
		## purpose of iterating through cavi B times

		cavi = q_gamma(X, y, sigma[i] sa[i], logodds[i], alpha[,i], mu[,i]) 

		#####       ----   store updates for i-th iteration    ----         ####

		logw[i]    = cavi$vlb       # variational lower bound

		## hyperparameter updates -- have to explicitly ask for update via input
		sigma[i]   = cavi$sigma     # changes iff update.sigma == TRUE
		sa[i]      = cavi$sa        # changes iff update.sa    == TRUE

		# variational parameter updates
		alpha[,i]  = cavi$alpha
		mu[,i]     = cavi$mu
		s[,i]      = cavi$s

		# mu.cov[,i] = cavi$post_mean # posterior mean for "additional" covars
		# note: for now, we leave out storing this variable -- doesn't affect
		# the calculations for posterior mean/variance that are described in the
		# simplest model		

	} # end of main for() -- logw, alpha, mu are updated, ready for PIP calc.


	# (normalized) weights, PIP, betas (all p-dim vectors)
	w    = normalizeLogWeights(logw)  #                    (B x 1)
	pip  = c(alpha %*% w)             # (p x B) (B x 1) -> (p x 1)
    beta = c((alpha * mu) %*% w)      # (n x B) (B x 1) -> (p x 1)

	# Prepare output
	output = list(pip = pip, beta = beta)
	return(output)

} # end of bvs() function







