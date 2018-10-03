

## calculate.R ----

## functions for common calculations


## function list:
	# l2()        : l2-norm of vector
	# var_ss()    : variance of var drawn from spike/slab prior (normal)
	# 


## specific to Carbonetto-Stephens Bayesian Variable Selection

	# varLB() : calculate the variational estimate to the marginal log-like
	# 



eps = .Machine$double.eps # machine precision



##### ---- l2() ---- #### 
# input: x vector
# output: 2-norm of the vector x (scalar)
l2 = function(x) {
	return(sqrt(sum(x * x)))
} # end of l2() function


##### ---- runif_mat() ---- #### 
# input:  m  : # of rows in matrix
#         n  : # of cols in matrix
# output: (m x n) matrix with each element ~ unif(0, 1)
runif_mat = function (m, n) {
	return(matrix(runif(m * n), m, n))
} # end of runif_mat() function




#### ---- var_ss() ---- ####
# variance of a random variable from from a spike and slab prior
# Var(X), where X ~ p N(x|mu, sigma) + (1 - p) delta0(x)
# X drawn from normal with probability p, X is 0 with probability 1 - p
# input: p     : probability that X is drawn from normal
#        mu    : mean parameter of the normal distribution
#        sigma : variance parameter of the normal distribution
# derive using law of total variance
var_ss = function(p, mu, sigma) {
	return(p * (sigma + mu^2) - (p * mu)^2)
} # end of var_ss() function


#### ---- logpexp ---- ####
# logpexp = log(1 + exp(x))
logpexp = function(x) {
	y = x 
	i = which(x < 16) # for large values of x, log(1 + exp(x)) = x (approx.)
	y[i] = log(1 + exp(x[i]))
	return(y)
}


#### ---- sigmoid() ---- ####
sigmoid = function(x) {
	return(1 / (1 + exp(-x)))
} # end of sigmoid() function


#### ---- logit() ---- ####
logit = function(x) {
	return(log((x + eps) / ((1 - x) + eps)))
} # end of sigmoid() function




#### ---- logsigmoid() ---- ####
# same as log(sigmoid(x))
logsigmoid = function() {
	return(-logpexp(-x))
} # end of logsigmoid() function




############## Carbonetto-Stephens Bayesian Variable Selection #################


##### ---- varLB() ---- #### 
# calculate the variational estimate to the marginal log-like
# analytical expression can be found in Carbonetto, Stephens (2012), page 81
# input: Xr    : X %*% (alpha * mu)                               (n x 1)
#        d     : is diagonal of X'X                               (n x 1)
#        y     : response vector                                  (n x 1)
#        sigma : VARIANCE of the RESIDUAL                         (1 x 1)
#        alpha : vector of prob that i-th coef ~ N(mu[i], s[i])   (p x 1)
#        mu    : MEAN VECTOR of the COEFFICIENTS                  (p x 1)
#        s     : VARIANCE VECTOR of the COEFFICIENTS              (p x 1)
#        sa    : sigma * sa is prior variance of coefficients     (1 x 1)
##
# note: logsigmoid(logodds) = log(p), where p is the prior probability of
#       beta_k ~ N(0, sigma * sa)
varLB = function(Xr, d, y, sigma, alpha, mu, s, logodds, sa) {f

	sigma_b = sigma * sa # prior variance of coefficients

	n = length(y)
	
	x0 = - n / 2 * log(2 * pi * sigma) -  l2(y - Xr)^2 / (2 * sigma) -
		   1 / (2 * sigma) * sum(d * var_ss(alpha, mu, s))
	
	x1 = sum((alpha - 1) * logodds) + logsigmoid(logodds)
	
	x2 = (sum(alpha) + sum(alpha * log(s / sigma_b)) - 
		      sum(alpha / sigma_b * (s + mu^2))) / 2 -
	      sum(alpha * log(alpha + eps)) - 
	      sum((1 - alpha) * log(1 - alpha + eps))


	return(x1 + x2 + x3)

} # end of varLB() function




############## Carbonetto-Stephens Bayesian Variable Selection #################




# end of calulate.R file
