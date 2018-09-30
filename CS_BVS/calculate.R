

## calculate.R ----

## functions for common calculations


## function list:
	# l2()    : l2-norm of vector
	# 
	# 


## specific to Carbonetto-Stephens Bayesian Variable Selection

	# varLB() : calculate the variational estimate to the marginal log-like
	# 



##### ---- l2() ---- #### 
# input: x vector
# output: 2-norm of the vector x (scalar)
l2 = function(x) {
	return(sqrt(sum(x * x)))
} # end of l2() function




##### ---- varLB() ---- #### 
# calculate the variational estimate to the marginal log-like
# input: Xr    : X %*% (alpha * mu)                             (n x 1)
#        d     : is diagonal of X'X t                           (n x 1)
#        y     : response vector                                (n x 1)
#        sigma : VARIANCE of the RESIDUAL                       (1 x 1)
#        alpha : vector of prob that i-th coef N(mu[i], s[i])   (p x 1)
#        mu    : MEAN VECTOR of the COEFFICIENTS                (p x 1)
#        s     : VARIANCE VECTOR of the COEFFICIENTS            (p x 1)
varLB = function(Xr, d, y, sigma, alpha, mu, s) {

	n = length(y)
	- n / 2 * log(2 * pi * sigma) -  

} # end of varLB() function
