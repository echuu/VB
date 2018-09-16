
## componentVB.R

## first algorithm described in Huang, Wang, Liang paper
## component-wise variational bayes

library(boot) # for inv.logit() function

n = 40
p = 8   # dimension of regression coefficients

## true parameters
beta0 = as.matrix(c(3, 1.5, 0, 0, 2, 0, 0, 0))
sigma = 3 

# specify hyperparameters: a0, b0, nu, lambda
# v1 in the spike & slab prior should also be pre-specified?


#######################    intialize parameters    #############################

# parameter vectors
mu_b         = numeric(p)
sigmasq_b    = numeric(p)
phi          = numeric(p)

# intial values for theta, sigmasq (both scalars)
theta   = 0.5
sigmasq = 1
#######################    intialize parameters    #############################


# init: (mu_1, ..., mu_p), (sigmasq_1, ..., sigmasq_p), (phi_1, ..., phi_p),
#       theta_hat, sigmasq_hat


# iteratively update the parameters (1 to p)
while (!converge) {
  # update beta here? 
  beta_bar = phi * mu_b
  for (j in 1:p) {
    
    ## update q_j (fully specified by mu_j, sigmasq_j)
    ## *** need to pre-specify v1 ***
    mu_b[j] = (t(y - X[,-j] %*% beta_bar[-j]) %*% X[,j]) / (n + 1 / v1)
    sigmasq_b[j] = sigmasq / (n + 1 / v1) 
    
    ## update phi_j
    phi[j] = inv.logit(logit(theta) - 0.5 * log(v1 * sigmasq / sigmasq_b[j]) + 
                         mu_b[j]^2 / (2 * sigmasq_b[j]))
    
  } # end for

  ## update theta_hat
  theta = (sum(phi) + a0 - 1) / (p + a0 + b0 - 2)
  
  ## update sigmasq
  num = sum((y - X %*% beta_bar)^2) + 
    sum((n * (1 - phi) + 1 / v1) * phi * mu_b^2 + 
          (n + 1 / v1) * phi * sigmasq_b) + nu * lambda
  sigmasq = num / (n + prod(phi) + nu + 2)
  
  # check for convergence
  
} # end while

# display parameter estimates


