
# errorMessages.R

# contains functions for checking validity of input throughout the algorithm

ARG_MISSING_WARNING   = "Missing arguments -- all arguments must be given."
MISSING_VALUE_WARNING = "Input should have numeric with no missing values"
DIMENSION_WARNING     = "Dimensions do not match"

checkInitialInput = function(X, y, sigma, sa, logodds, n, p) {
	# check that X, y, sigma, sa, logodds are all provided
	ARG_MISSING = missing(X)  || missing (y) || missing(sigma) || 
                  missing(sa) || missing(logodds)
    if (ARG_MISSING) {
    	stop(ARG_MISSING_WARNING)
    }


    # check for missing values in X, y
    # check that dimension of y matches the # of rows of X
    X_ERR     = !(is.matrix(X) & is.numeric(X) & sum(is.na(X)) == 0)
	Y_ERR     = !is.numeric(y) | sum(is.na(y)) > 0
	Y_DIM_ERR = length(y) != n

	if (X_ERR | Y_ERR) {
		stop(MISSING_VALUE_WARNING)
	}
	if (Y_DIM_ERR) {
		stop(DIMENSION_WARNING)
	}
	
} # end of checkInitialInput() function



# check that after processing the hyperparameters, they are:
# logodds: 1 x B or p x B
# sigma  : B x 1
# sa     : B x 1
checkHyperDims = function(logodds, sigma, sa, B) {
	if (length(sigma) != B | length(sa) != B | ncol(logodds) != B) {
		stop("Dimensions of hyperparameters is incompatible")
	} 
} # end of checkHyperDims() function







