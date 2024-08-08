#sd <- function(K, r,r1, n) {
#P0 <- P0 + (((1 - r1^K) / (1 - r1)) + ((r * r1^(K - 1)) / (1 - r)))
  #ps <- take any positive value for ps...
  Pn <- rep(0, ps)
  N <- rep(0, ps)
  for (i in 0:ps) {
    if (i < K) {
      P <- (r1^i) * P0
    } else {
      P <- ((r1^(K - 1)) * (r^(i - K + 1))) * P0
    }
    #c<- take any initial value for c to store value in matrix form
    Pn[c] <- P
    N[c] <- i
    c <- c + 1
  }
  
  #T<- take any positive value to generate the random sample
  iter <- 5000 #increase the value of iteration to get more accuracy in the estimated value
  # for mle use the optimization method 
  # for Bayesian estimates use different priors
 # Pr_r <- rbeta(T,1,3) #change the prior according to the condition of the parameter
#  Pr_r1 <- rgamma(T,1,3)
  while (init < n) {
    #R <- sample(N, T, replace = TRUE, prob = Pn) #for random sample generation change the function according to the conditions of density 
    #V_r <- rep(0, iter)
    #V_r1 <- rep(0, iter)
    for (j in 1:T) {
      r_j <- Pr_r[j]
      for (r in 1:iter) {
        r1_r <- Pr_r1[r]
        #P0 <- P0 + (((1 - r1_r^K) / (1 - r1_r)) + ((r_j * r1_r^(K - 1)) / (1 - r_j)))
        
        E1 <- 1
        for (i in 1:T) {
          i1 <- R[i]
          if (i1 < K) {
            E1a <- (r1_r^i1) * P0
            E1 <- E1 * E1a
          } else {
            E2a <- ((r1_r^(K - 1)) * (r_j^(i1 - K + 1))) * P0
            E1 <- E1 * E2a
          }
        }
        V_r[j] <- E1
        V_r1[r] <- E1
      }
    }
    
    Po_r <- sample(Pr_r, Tasir, replace = TRUE, prob = V_r)
    Po_r1 <- sample(Pr_r1, Tasir, replace = TRUE, prob = V_r1)
    
    Pr_r <- Po_r
    Pr_r1 <- Po_r1
  }
  
  #for an estimate of r
  est <- mean(Po_r)
   #for r1
  est1 <- mean(Po_r1)
  # For the the value of square error, risk, confidence, and credible interval generate the formula accordingly for parameters
  
  
}

sd(3,0.3,0.5,50)
