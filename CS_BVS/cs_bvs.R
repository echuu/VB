#
# Compute fully-factorized variational approximation for Bayesian
# variable selection in linear (family = "gaussian") or logistic
# regression (family = "binomial"). See varbvs.Rd for details.
varbvs <- function (X, Z, y, family = c("gaussian","binomial"), sigma, sa,
                    logodds, weights, resid.vcov, alpha, mu, eta, 
                    update.sigma, update.sa, optimize.eta, initialize.params,
                    update.order, nr = 100, sa0 = 1, n0 = 10, tol = 1e-4,
                    maxiter = 1e4, verbose = TRUE) {

    # Get the number of samples (n) and variables (p).
    n <- nrow(X)
    p <- ncol(X)

    # (1) CHECK INPUTS
    # ----------------
    
    # Get candidate settings for the variance of the residual (sigma),
    # if provided. Note that this option is not valid for a binary trait.
    if (missing(sigma)) {
      sigma <- var(y)
      update.sigma.default <- TRUE
    } else {
      sigma <- c(sigma)
      update.sigma.default <- FALSE
      if (family == "binomial")
        stop("Input sigma is not allowed for family = binomial")
    }
    
    # Get candidate settings for the prior variance of the coefficients,
    # if provided.
    if (missing(sa)) {
      sa <- 1
      update.sa.default <- TRUE
    } else {
      sa <- c(sa)
      update.sa.default <- FALSE
    }

    # Get candidate settings for the prior log-odds of inclusion. This
    # may either be specified as a vector, in which case this is the
    # prior applied uniformly to all variables, or it is a p x ns
    # matrix, where p is the number of variables and ns is the number of
    # candidate hyperparameter settings, in which case the prior
    # log-odds is specified separately for each variable. A default
    # setting is only available if the number of other hyperparameter
    # settings is 1, in which case we select 20 candidate settings for
    # the prior log-odds, evenly spaced between log10(1/p) and -1.
    if (missing(logodds)) {
      if (length(sigma) == 1 & length(sa) == 1)
        logodds <- seq(-log10(p),-1,length.out = 20)
      else
        stop("logodds can only be missing when length(sigma) = length(sa) = 1")
    }
    if (!is.matrix(logodds)) {
      prior.same <- TRUE
      logodds    <- t(matrix(logodds))
    } else if (nrow(logodds) != p) {
      prior.same <- TRUE
      logodds    <- t(matrix(logodds))
    } else
      prior.same <- FALSE

    # Here is where I ensure that the numbers of candidate hyperparameter
    # settings are consistent.

    # ns is the number of times the main for loop iterates
    # = the number of hyperparameters
    # this number dictates the number of columns in the storage data structures
    ns <- max(length(sigma),length(sa),ncol(logodds))
    
    if (length(sigma) == 1)
      sigma <- rep(sigma,ns)
    if (length(sa) == 1)
      sa <- rep(sa,ns)
    if (ncol(logodds) == 1)
      logodds <- rep.col(logodds,ns)
    if (length(sigma) != ns | length(sa) != ns | ncol(logodds) != ns)
      stop("options.sigma, options.sa and options.logodds are inconsistent")

    # Set initial estimates of variational parameter alpha.
    initialize.params.default <- TRUE
    if (missing(alpha)) {
      alpha <- rand(p,ns)
      alpha <- alpha / rep.row(colSums(alpha),p)
    } else
      initialize.params.default <- FALSE
    if (!is.matrix(alpha))
      alpha <- matrix(alpha)
    if (nrow(alpha) != p)
      stop("Input alpha must have as many rows as X has columns")
    if (ncol(alpha) == 1)
      alpha <- rep.col(alpha,ns)

    # Set initial estimates of variational parameter mu.
    if (missing(mu))
      mu <- randn(p,ns)
    else
      initialize.params.default <- FALSE    
    if (!is.matrix(mu))
      mu <- matrix(mu)
    if (nrow(mu) != p)
      stop("Input mu must have as many rows as X has columns")
    if (ncol(mu) == 1)
      mu <- rep.col(mu,ns)

    # Determine whether to find a good initialization for the
    # variational parameters.
    if (missing(initialize.params))
      initialize.params <- initialize.params.default
    else if (initialize.params & ns == 1)
      stop(paste("initialize.params = TRUE has no effect when there is",
                 "only one hyperparameter setting"))


    # (3) PREPROCESSING STEPS
    # -----------------------
    if (family == "gaussian") {
      
        # Adjust the inputs X and y so that the linear effects of the
        # covariates (Z) are removed. This is equivalent to integrating
        # out the regression coefficients corresponding to the covariates
        # with respect to an improper, uniform prior.
        out <- remove.covariate.effects(X,Z,y)
        X   <- out$X
        y   <- out$y
        rm(out)

    } # end of if (gaussian)

    # Add row and column names to X if they are not provided.
    if (is.null(rownames(X)))
        rownames(X) <- 1:n
    if (is.null(colnames(X)))
        colnames(X) <- paste0("X",1:p)
    
    # (4) INITIALIZE STORAGE FOR THE OUTPUTS
    # --------------------------------------
    # Initialize storage for the variational estimate of the marginal
    # log-likelihood for each hyperparameter setting (logw), and the
    # variances of the regression coefficients (s), and the posterior
    # mean estimates of the coefficients for the covariates (mu.cov),
    # which includes the intercept.
    logw   <- rep(0,ns)
    s      <- matrix(0,p,ns)
                     
    # (5) FIT BAYESIAN VARIABLE SELECTION MODEL TO DATA
    # -------------------------------------------------
    if (ns == 1) {

        # Find a set of parameters that locally minimize the K-L
        # divergence between the approximating distribution and the exact
        # posterior.

        out   <- outerloop(X,Z,y,family,weights,resid.vcov,SZy,SZX,c(sigma),
                           c(sa),c(logodds),c(alpha),c(mu),c(eta),update.order,
                           tol,maxiter,verbose,NULL,update.sigma,update.sa,
                           optimize.eta,n0,sa0)
        logw     <- out$logw
        sigma    <- out$sigma
        sa       <- out$sa
        alpha[]  <- out$alpha
        mu[]     <- out$mu
        s[]      <- out$s
        eta[]    <- out$eta

    } else {

        # If a good initialization isn't already provided, find a good
        # initialization for the variational parameters. Repeat for each
        # candidate setting of the hyperparameters
        if (initialize.params) {

            # Repeat for each setting of the hyperparameters.
            # the following for loop appears to do the same thing as the for
            # loop that comes directly after this if-statement, but:
            # the assingnments that come right after this for-loop take the 
            # choice of variational parameters (from the pool of 1:ns params)
            # that give the largest variational lower bound and uses these 
            # for ALL initial values for the corresponding variational parameter

            for (i in 1:ns) {
                out <- outerloop(X,Z,y,family,weights,resid.vcov,SZy,SZX,sigma[i],
                                 sa[i],logodds[,i],alpha[,i],mu[,i],eta[,i],
                                 update.order,tol,maxiter,verbose,i,update.sigma,
                                 update.sa,optimize.eta,n0,sa0)
                
                logw[i]    <- out$logw

                sigma[i]   <- out$sigma
                sa[i]      <- out$sa

                alpha[,i]  <- out$alpha
                mu[,i]     <- out$mu
                s[,i]      <- out$s

            } # end for()

            # Choose an initialization common to all the runs of the
            # coordinate ascent algorithm. This is chosen from the
            # hyperparameters with the highest variational estimate of the
            # marginal likelihood.
            i     <- which.max(logw)        # find iteration w/ largest VLB
            alpha <- rep.col(alpha[,i],ns)  # repeat for ALL initial alphas
            mu    <- rep.col(mu[,i],ns)     # repeat for ALL initial mus
            
            if (update.sigma)               # hyperparameters update
              sigma = rep(sigma[i],ns)      # update variance of residual
                                            # all hyperparams same (??) why
            if (update.sa)                 
              sa = rep(sa[i],ns)            # update (scaled) variance of coeffs


              ## if both of the if() statements above execute, then every iter
              ## of the main loop below will use the same variational parameters
              ## and the same hyperparameters, which doesn't quite make sense


        } # end of if (initialize.param) chunk

        # Compute a variational approximation to the posterior distribution
        # for each candidate setting of the hyperparameters.
        
        # For each setting of the hyperparameters, find a set of
        # parameters that locally minimize the K-L divergence between the
        # approximating distribution and the exact posterior.
        for (i in 1:ns) {
            out <- outerloop(X, Z, y, family, weights, resid.vcov, SZy, SZX, 
                             sigma[i],sa[i], logodds[,i], alpha[,i], mu[,i], 
                             update.order, tol, maxiter, verbose, i,
                             update.sigma, update.sa, optimize.eta, n0, sa0)
            logw[i]    <- out$logw
            sigma[i]   <- out$sigma
            sa[i]      <- out$sa


            # for each iteration, we update the hyperparameters
            # phi_i = (alpha_i, mu_i, s_i) using CAVI, which happens when we
            # call outerloop()
            # previous for loop that performs the 'ns' iterations is to find
            # good starting values phi_i so that this for loop has faster
            # convergence (but doesn't the previous for loop do the work of 
            # this loop?) --- local convergence, so it might be better to start
            # with better initial values to get faster convergence

            alpha[,i]  <- out$alpha
            mu[,i]     <- out$mu
            s[,i]      <- out$s
        }

    } # end of if-else() chunk

    # (6) CREATE FINAL OUTPUT
    # -----------------------
    # Compute the normalized importance weights and the posterior
    # inclusion probabilities (PIPs) and mean regression coefficients
    # averaged over the hyperparameter settings.
    if (ns == 1) {
        w        <- 1
        pip      <- c(alpha)
        beta     <- c(alpha*mu)
        beta.cov <- c(mu.cov)
    } else {
        w        <- normalizelogweights(logw)
        pip      <- c(alpha %*% w)
        beta     <- c((alpha*mu) %*% w)
        beta.cov <- c(mu.cov %*% w)
    }

    if (prior.same)
        fit$logodds <- c(fit$logodds)
    else
        rownames(fit$logodds) <- colnames(X)

    
    return(fit)
}


