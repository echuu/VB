
# bvs.R 

## Bayesian Variable Selection via the method described in Carbonetto & Stephens
## Each iteration gives an alpha, which are used to generate an MCMC estimate of
## the posterior inclusion probability of the regression coefficients


## hyperparameters: sigma, sa
## variational parameters: alpha, mu, s


## input:
#          X       : (n x p) design matrix, n is # of obs, p is # of variables
#          sigma   : (1 x 1) variance of the RESIDUAL
#          sa      : (1 x 1) PRIOR variance of regression coefficients
#          logodds : (p x 1) prior log-odds of inclusion for each variable
#          n_ps    : (1 x 1) number of samples from the posterior 

## output:
#          wts       : (B x 1) normalized weights
#          pip       : (p x 1) probability of inclusion for each of the p coeffs
#          beta      : (p x 1) expected value of each of the betas



library(dplyr)


bvs = function(X, y, sigma, sa, logodds,
	           update.sigma, update.sa, n_ps = 100, sa0 = 1, n0 = 10,
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

    # resulting X will just be the centered version, and y = y - mean(y)
    new_Xy = processXy(X, y)
    X      = new_Xy$X
    y      = new_Xy$y 

    rm(new_Xy)

	#### ----          preprocess output and hyperparameters           ---- ####
	# --------------------------------------------------------------------------

	# B = # of times the main for loop iterates
    #   = the number of CANDIDATE hyperparameters
    # this number dictates the number of columns in the storage data structures
    # (defined as ns in C&S code)

	B = max(length(sigma), length(sa), ncol(logodds))

	hyper   = processHyper(logodds, sigma, sa, B, p)
	logodds = hyper$logodds            # 1 x B if input was scalar, else: p x B
	sigma   = hyper$sigma			   # B x 1 vector
	sa      = hyper$sa                 # B x 1 vector


	#### ----            initialize variational parameters             ---- ####
	# --------------------------------------------------------------------------
	alpha = runif_mat(p, B)    # (p x B) estimates of each iter stored col-wise
	mu    = runif_mat(p, B)    # (p x B) estimates of each iter stored col-wise	
	logw  = rep(0, B)                   # var. estimate of marginal log-like
	s     = matrix(0, p, B)             # variance of reg coeffs, stored col-wise

	

    #### ----                       outer loop                         ---- ####
	# --------------------------------------------------------------------------
	
	# (1) If initialization of the parameters is not provided in the input, 
	#     we iterate through the variational updates to find a better 
	#     initialization for the variational parameters. The values we use to
	#     run CAVI for this chunk are generated randomly when we defined above
	for (i in 1:B) {


		## first run through of algorithm to find optimal variational
		## parameter settings for more accurate posterior estimates
		cavi_0 = q_gamma(X, y, sigma[i], sa[i], logodds[,i], alpha[,i], mu[,i])

		logw[i]   = cavi_0$logw

		sigma[i]  = cavi_0$sigma
		sa[i]     = cavi_0$sa

		alpha[,i] = cavi_0$alpha
		mu[,i]    = cavi_0$mu
		s[,i]     = cavi_0$s

		## find the iteration with parameters that maximized the VLB
		## store these parameters in alpha, mu, sigma (*), sa (*)
		## (*): if update.sigma == TRUE
		i_max = which.max(logw)
		alpha = rep.col(alpha[, i_max], B) # repeat optimal alpha for all inits
		mu    = rep.col(mu[, i_max], B)    # repeat optimal mu for all inits


		if (update.sigma) {
			sigma = rep(sigma[i], B) # repeat same sigma for all B iterations
		}                            # doesn't result in identitcal results for
		                             # each of the B iterations in the
		                             # subsequent for loop ??

		if (update.sa) {
			sa = rep(sa[i], B)       # same concern/question as above
		}                             

	} # end of initilization for()


	# (2) Using the initialization chosen in step (1), we run CAVI.
	for (i in 1:B) {    # beginning of the outer loop of algorithm

		# inner loop -- optimize var lower bound

		## 10/5 questions/confusions:
		## details of these updates still need to be figured out..
		## seems that if we have initialize.params == TRUE, then all starting
		## values of alpha[,i], mu[,i] will all be the same, which would result
		## in the same updates stored in each iteration, which would defeat the
		## purpose of iterating through cavi B times

		## 10/8 update: 
		## we do want to be using the same INITIALIZATION for the variational
		## parameters, but these updates will change because the 
		## hyperparameter settings are different for each of the B iterations
		## (see 2nd paragraph on p. 82 of C&S (2012))
		## since they are sampled from an importance density (where this done?)

		## follow up question:
		## where in this process have we sampled the hyperparmeters from the
		## importance density mentioned in the paper? 

		## 10/9 update:
		## paper says that we use grid points for candidate values for theta
		## since the values of the variational parameters depend on these values
		## of theta, using the same initialization for the variational params
		## (which is done in the initial run), will not necessarily result
		## in the exact same cavi outputs, due to convergence to local optima


		cavi = q_gamma(X, y, sigma[i], sa[i], logodds[i], alpha[,i], mu[,i]) 

		#####       ----   store updates for i-th iteration    ----         ####

		logw[i]    = cavi$vlb       # variational lower bound

		## hyperparameter updates -- have to explicitly ask for update via input
		sigma[i]   = cavi$sigma     # changes iff update.sigma == TRUE
		sa[i]      = cavi$sa        # changes iff update.sa    == TRUE

		# variational parameter updates
		alpha[,i]  = cavi$alpha
		mu[,i]     = cavi$mu
		s[,i]      = cavi$s

	} # end of main for() -- logw, alpha, mu are updated, ready for PIP calc.


	# (normalized) weights, PIP, betas (all p-dim vectors)
	w    = normalizeLogWeights(logw)  #                    (B x 1)
	pip  = c(alpha %*% w)             # (p x B) (B x 1) -> (p x 1)
    beta = c((alpha * mu) %*% w)      # (n x B) (B x 1) -> (p x 1)

	# Prepare output
	output = list(wts = w, pip = pip, beta = beta)
	return(output)

} # end of bvs() function







