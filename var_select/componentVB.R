
## componentVB.R
## first algorithm described in Huang, Wang, Liang paper
## component-wise variational bayes

n = 40
p = 8   # dimension of regression coefficients

## true parameters
beta0 = as.matrix(c(3, 1.5, 0, 0, 2, 0, 0, 0))
sigma = 3 

# init: (mu_1, ..., mu_p), (sigmasq_1, ..., sigmasq_p), (phi_1, ..., phi_p),
#       theta_hat, sigmasq_hat

# iteratively update the parameters (1 to p)
while (!converge) {
  for (j in 1:p) {
    
    # update q_j (fully specified by mu_j, sigmasq_j)
    
    # update phi_j
    
  } # end for

  # update theta_hat
  
  # update sigma_sq
  
  # check for convergence
  
} # end while

# display parameter estimates


