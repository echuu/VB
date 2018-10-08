
Implementation of various Bayesian Variable Selection Techniques

CB_BVS: implementation of BVS method discussed in Carbonetto and Stephens (2012)

Descriptions of files/scripts:

        bvs.R           : the main routine of the algorithm, takes user input, 
                          the primary purpose is to run the outer loop of the
                          algorithm, making calls to q_gamma() every iteration
                          to retrieve the cavi updates for the variational
                          params; the posterior inclusion probability, as well
                          as the expected value of each of the coefficients are
                          returned at the end of this routine

        q_gamma.R       : inner loop of the variational bvs algorithm, performs
                          the cavi updates for each of the variational params
                          until convergence, returns correspoding estimates 

        vbUpdate.R      : performs ONE iteration of the CAVI update; this func.
                          is called within q_gamma()'s  iterating scheme                                                


        processXy.R     : centers X matrix, subtracts mean of response from y

        calculate.R     : contains various calculations/procedures that are used
                          throughout the algorithm
      
        errorMessages.R : stores and prints the error messages corresponding to
                          various checks throughout the algorithm


        cs_bvs.R        : C&S implementation of the main variational bvs alg;
                          returns the posterior inclusions probs of the coeffs

        cs_outerloop.R  : performs the outer loop described in the alg. outline
                          in C&S (2012). In my implementation, this step is
                          written directly within the outer loop.



