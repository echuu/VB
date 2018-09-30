    

# vbUpdate() performs ONE iteration of the CAVI update to
# maximize the variational lower bound for Bayesian variable 
# selection in linear regression.

## input:
#          X       : (n x p) design matrix, n is # of obs, p is # of variables
#          sigma   : (1 x 1) variance of the residual
#          sa      : (1 x 1) prior variance of regression coefficients
#          logodds : (p x 1) prior log-odds of inclusion for each variable 
#          xy      : (p x 1) X'y
#          d       : (p x 1) diagonal entires of X'X matrix
#          alpha0  : (p x 1) curr param of var approx of PIP(i) = alpha0[i]
#          mu0     : (p x 1) mean of coeff_i given that it's in model, mu0[i]
#          Xr0     : (n x 1) X * (alpha0 * mu0)

## output: updated variational parameters
#          alpha   : (p x 1) curr param of var approx of PIP(i) = alpha[i]
#          mu      : (p x 1) mean of coeff_i given that it's in model, mu[i]
#          Xr      : (n x 1) X * (alpha * mu)


vbUpdate = function(X, sigma, sa, logodds, xy, d, alpha0, mu0, Xr0) {


    # order of variational updates -- default: 1:p
    updates = 1:length(alpha)

    # initialize parameters
    alpha = alpha0
    mu    = mu0
    Xr    = Xr0

    # CAVI updates 
    for (j in updates) {

        # posterior variance
        s = sa * sigma / (sa * d[j] + 1)

        # posterior mean.
        r     = alpha[j] * mu[j]
        mu[j] = s / sigma * (xy[j] + d[j] * r - sum(X[, j] * Xr))

        # variational estimate of the posterior inclusion probability.
        logBF    = (log(s / (sa * sigma)) + mu[j]^2 / s) / 2 # update mistake ?
        alpha[j] = sigmoid(logodds[j] + logBF)

        # Update Xr = X*r.
        Xr = Xr + (alpha[j] * mu[j] - r) * X[,j]
        
    } # end of update loop


    return(list(alpha = alpha, mu = mu, Xr = Xr))

} # end of vbUpdate() function
