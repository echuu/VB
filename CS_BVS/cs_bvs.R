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
  # Check input Z.
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z) | sum(is.na(Z)) > 0)
      stop("Input Z must be a numeric matrix with no missing values.")
    if (nrow(Z) != n)
      stop("Inputs X and Z do not match.")
    storage.mode(Z) <- "double"
  }
  
  # Add intercept.
  if (is.null(Z))
    Z <- matrix(1,n,1)
  else
    Z <- cbind(1,Z)

  # Get choice of regression model.
  family <- match.arg(family)

  # (2) PROCESS OPTIONS
  # -------------------
  if (!is.finite(maxiter))
    stop("Input maxiter must be a finite number")
  
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
  ns <- max(length(sigma),length(sa),ncol(logodds))
  if (length(sigma) == 1)
    sigma <- rep(sigma,ns)
  if (length(sa) == 1)
    sa <- rep(sa,ns)
  if (ncol(logodds) == 1)
    logodds <- rep.col(logodds,ns)
  if (length(sigma) != ns | length(sa) != ns | ncol(logodds) != ns)
    stop("options.sigma, options.sa and options.logodds are inconsistent")


  # Set the weights---used to indicate that different observations
  # have different variances---or the covariance matrix of the
  # residual (resid.vcov). Note that only one of the weights and
  # residual covariance matrix can be non-NULL.
  if (missing(weights))
    weights <- NULL
  if (missing(resid.vcov)) 
    resid.vcov <- NULL
  

  ### **** NOTE: if resid.vcov not passed in, then it is SET TO NULL HERE ****##


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
    SZy <- out$SZy
    SZX <- out$SZX
    rm(out)

  } # end of if (gaussian)

  # Add row and column names to X if they are not provided.
  if (is.null(rownames(X)))
    rownames(X) <- 1:n
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X",1:p)

  # Add column names to Z if they are not already provided.
  if (is.null(colnames(Z)) & ncol(Z) > 1)
    colnames(Z) <- c("(Intercept)",paste0("Z",1:(ncol(Z) - 1)))
  else
    colnames(Z)[1] <- "(Intercept)"
  
  # (4) INITIALIZE STORAGE FOR THE OUTPUTS
  # --------------------------------------
  # Initialize storage for the variational estimate of the marginal
  # log-likelihood for each hyperparameter setting (logw), and the
  # variances of the regression coefficients (s), and the posterior
  # mean estimates of the coefficients for the covariates (mu.cov),
  # which includes the intercept.
  logw   <- rep(0,ns)
  s      <- matrix(0,p,ns)
  mu.cov <- matrix(0,ncol(Z),ns)
                   
  # (5) FIT BAYESIAN VARIABLE SELECTION MODEL TO DATA
  # -------------------------------------------------
  if (ns == 1) {

    # Find a set of parameters that locally minimize the K-L
    # divergence between the approximating distribution and the exact
    # posterior.

    out      <- outerloop(X,Z,y,family,weights,resid.vcov,SZy,SZX,c(sigma),
                          c(sa),c(logodds),c(alpha),c(mu),c(eta),update.order,
                          tol,maxiter,verbose,NULL,update.sigma,update.sa,
                          optimize.eta,n0,sa0)
    logw     <- out$logw
    sigma    <- out$sigma
    sa       <- out$sa
    mu.cov[] <- out$mu.cov
    alpha[]  <- out$alpha
    mu[]     <- out$mu
    s[]      <- out$s
    eta[]    <- out$eta
    if (verbose)
      cat("\n")
  } else {

    # If a good initialization isn't already provided, find a good
    # initialization for the variational parameters. Repeat for each
    # candidate setting of the hyperparameters.
    if (initialize.params) {
      if (verbose) {
        cat("Finding best initialization for",ns,"combinations of",
            "hyperparameters.\n");
        cat("-iteration-   variational    max.   incl variance params\n");
        cat("outer inner   lower bound  change   vars   sigma      sa\n");
      }

      # Repeat for each setting of the hyperparameters.
      for (i in 1:ns) {
        out <- outerloop(X,Z,y,family,weights,resid.vcov,SZy,SZX,sigma[i],
                         sa[i],logodds[,i],alpha[,i],mu[,i],eta[,i],
                         update.order,tol,maxiter,verbose,i,update.sigma,
                         update.sa,optimize.eta,n0,sa0)
        logw[i]    <- out$logw
        sigma[i]   <- out$sigma
        sa[i]      <- out$sa
        mu.cov[,i] <- out$mu.cov
        alpha[,i]  <- out$alpha
        mu[,i]     <- out$mu
        s[,i]      <- out$s
        eta[,i]    <- out$eta
      }
      if (verbose)
        cat("\n")

      # Choose an initialization common to all the runs of the
      # coordinate ascent algorithm. This is chosen from the
      # hyperparameters with the highest variational estimate of the
      # marginal likelihood.
      i     <- which.max(logw)
      alpha <- rep.col(alpha[,i],ns)
      mu    <- rep.col(mu[,i],ns)
      if (optimize.eta)
        eta <- rep.col(eta[,i],ns)
      if (update.sigma)
        sigma <- rep(sigma[i],ns)
      if (update.sa)
        sa <- rep(sa[i],ns)
    }

    # Compute a variational approximation to the posterior distribution
    # for each candidate setting of the hyperparameters.
    
    # For each setting of the hyperparameters, find a set of
    # parameters that locally minimize the K-L divergence between the
    # approximating distribution and the exact posterior.
    for (i in 1:ns) {
      out <- outerloop(X,Z,y,family,weights,resid.vcov,SZy,SZX,sigma[i],sa[i],
                       logodds[,i],alpha[,i],mu[,i],eta[,i],update.order,tol,
                       maxiter,verbose,i,update.sigma,update.sa,optimize.eta,
                       n0,sa0)
      logw[i]    <- out$logw
      sigma[i]   <- out$sigma
      sa[i]      <- out$sa
      mu.cov[,i] <- out$mu.cov
      alpha[,i]  <- out$alpha
      mu[,i]     <- out$mu
      s[,i]      <- out$s
      eta[,i]    <- out$eta
    }
    if (verbose)
      cat("\n")
  }

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

  if (family == "gaussian") {
    fit <- list(family = family,n0 = n0,sa0 = sa0,mu.cov = mu.cov,
                update.sigma = update.sigma,update.sa = update.sa,
                prior.same = prior.same,optimize.eta = FALSE,logw = logw,
                w = w,sigma = sigma,sa = sa,logodds = logodds,alpha = alpha,
                mu = mu,s = s,eta = NULL,pip = pip,beta = beta,
                beta.cov = beta.cov,y = y)
    class(fit) <- c("varbvs","list")

    if (is.null(weights) & is.null(resid.vcov) & ncol(Z) == 1) {
        
      # Compute the proportion of variance in Y---only in the
      # unweighted, i.i.d. case when there are no additional covariates
      # included in the model.
      if (verbose)
        cat("Estimating proportion of variance in Y explained by model.\n");
      fit$model.pve <- varbvspve(fit,X,nr)

      # Compute the proportion of variance in Y explained by each
      # variable. This can only be estimated in the i.i.d. case when
      # there are no additional covariates included in the model.
      fit$pve           <- matrix(0,p,ns)
      rownames(fit$pve) <- colnames(X)
      sx                <- var1.cols(X)
      for (i in 1:ns) 
        fit$pve[,i] <- sx*(mu[,i]^2 + s[,i])/var1(y)
    } else {
      fit$model.pve <- NULL
      fit$pve       <- NULL
    }

    # Restore the inputted X and y.
    X <- X + Z %*% SZX
    y <- y + c(Z %*% SZy)
  
    # Compute the fitted values for each hyperparameter setting.
    fit$fitted.values <- varbvs.linear.predictors(X,Z,family,mu.cov,alpha,mu)

    # Compute the residuals for each hyperparameter setting.
    fit$residuals <- y - fit$fitted.values
    fit$residuals.response
  } else if (family == "binomial") {
    fit <- list(family = family,n0 = n0,mu.cov = mu.cov,sa0 = sa0,
                update.sigma = FALSE,update.sa = update.sa,
                optimize.eta = optimize.eta,prior.same = prior.same,
                logw = logw,w = w,sigma = NULL,sa = sa,logodds = logodds,
                alpha = alpha,mu = mu,s = s,eta = eta,pip = pip,beta = beta,
                beta.cov = beta.cov,pve = NULL,model.pve = NULL)
    class(fit) <- c("varbvs","list")

    # Compute the fitted values for each hyperparameter setting.
    fit$fitted.values <-
      sigmoid(varbvs.linear.predictors(X,Z,family,mu.cov,alpha,mu))
    
    # Compute the "deviance" and "response" residuals for
    # hyperparameter setting.
    fit$residuals <-
      list(deviance = resid.dev.logistic(matrix(y,n,ns),fit$fitted.values),
           response = y - fit$fitted.values)
  }
  
  # Add names to some of the outputs.
  hyper.labels                <- paste("theta",1:ns,sep = "_")
  rownames(fit$alpha)         <- colnames(X)
  rownames(fit$mu)            <- colnames(X)
  rownames(fit$s)             <- colnames(X)
  names(fit$beta)             <- colnames(X)
  names(fit$pip)              <- colnames(X)
  rownames(fit$mu.cov)        <- colnames(Z)
  names(fit$beta.cov)         <- colnames(Z)
  rownames(fit$fitted.values) <- rownames(X)
  colnames(fit$mu.cov)        <- hyper.labels
  colnames(fit$alpha)         <- hyper.labels
  colnames(fit$mu)            <- hyper.labels
  colnames(fit$s)             <- hyper.labels
  colnames(fit$fitted.values) <- hyper.labels
  if (family == "gaussian") {
    rownames(fit$residuals) <- rownames(X)
    colnames(fit$residuals) <- hyper.labels
  } else {
    rownames(fit$eta)                <- rownames(X)
    colnames(fit$eta)                <- hyper.labels
    rownames(fit$residuals$deviance) <- rownames(X)
    rownames(fit$residuals$response) <- rownames(X)
    colnames(fit$residuals$deviance) <- hyper.labels
    colnames(fit$residuals$response) <- hyper.labels
  }
  if (prior.same)
    fit$logodds <- c(fit$logodds)
  else
    rownames(fit$logodds) <- colnames(X)
  if (!is.null(fit$pve))
    colnames(fit$pve) <- hyper.labels
  return(fit)
}

