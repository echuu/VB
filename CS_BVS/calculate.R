

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


# Replicate vector x to create an m x n matrix, where m = length(x).
##### ---- rep_col() ---- #### 
# input:  x  : m-dimensional vector
#         n  : # of times to replicate the column vector x
# output: (m x n) matrix where x is repeated column-wise n times
rep_col <- function (x, n) {
	return(matrix(x, length(x), n, byrow = FALSE))
} # end of rep_col() function


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


## input:  user input versions of logodds, sigma, sa;
##         B is the number of candidate hyperparameters, dictates dimensions
##         p is the dimension of the coefficient vector
## output: logodds (p x B) or (1 x B), sigma (1 x B), sa (1 x B)
processHyper = function(logodds, sigma, sa, B, p) {
	

	if (length(sigma) == 1) {
		sigma <- rep(sigma, B)
	}
	if (length(sa) == 1) {
		sa = rep(sa, B)
	}

	## all of the following if() statements will execute if logodds is a scalar
	## keep them separate for now, since we are (perhaps incorrectly?)
	## assuming that the vector that is passed in must be 1 x 1 if we are to
	## apply the prior uniformy across all the variables, i.e., if logodds is
	## not passed in as a matrix, then it must be a 1 x 1 vector (scalar)

	if (length(logodds) == 1) {
		logodds = rep(logodds, p)       
	}

	if (!ismatrix(logodds)) {           # if input for logodds is SCALAR
		logodds = t(matrix(logodds))    # prior is applied unif to all variables
	}                                   # use matrix form for syntax conformity

	if (ncol(logodds) == 1) {           # if input for logodds is SCALAR, then
		logodds = rep_col(logodds, B)   # prev. step will result in 1x1 matrix
	} 									# * logodds is (1 x B) after this step *


	checkHyperDims(logodds, sigma, sa, B, p) # check dimensions of hyperparams
	
	return(list(logodds, sigma, sa))
}


#### ---- normalizeLogWeights() ---- ####
# input:  logw : vector of log probabilities
# output: vector of normalized probabilities (sum to 1)
normalizeLogWeights = function(logw) {

	c = max(logw)
	w = exp(logw - c) # same as dividing by max on regular scale

	return(w / sum(w))
} # end of normalizeLogWeights()



############## Carbonetto-Stephens Bayesian Variable Selection #################




# end of calulate.R file
