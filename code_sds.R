#Program_Table1:
library(MLmetrics)  #for MSE  
 # Function to generate random sample from the given density
  s <- function(K, rho, rho1, n) {
    pr <- rep(0, n)
    for (i in 1:n) {
      if (i < K) {
        pr[i] <- (rho1^(i)) * (((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))^(-1)
      } else {
        pr[i] <- (rho1^(K - 1)) * (rho^(i - K + 1)) * ((((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))^(-1))
      }
    }
    sa <- sample(0:(n - 1), size = n, replace = TRUE, prob = pr)
    return(sa)
  }
  
 
  ll <- function(parameters, K, n, sa) {
    rho <- parameters[1]
    rho1 <- parameters[2]
    
    pr <- rep(0, n)
    for (i in 1:n) {
      if (sa[i] < K) {
        pr[i] <- (rho1^(sa[i])) * (((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))^(-1)
      } else {
        pr[i] <- (rho1^(K - 1)) * (rho^(sa[i] - K + 1)) * ((((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))^(-1))
      }
    }
    
    lld <- sum(log(pr))
    return(-lld)  
  }
  
  # To fix the output, you can use set.seed followed by any random value, for example: set.seed(123) or set.seed(134567). 
  true_rho <- 0.7
  true_rho1<- 0.5
  rho <- 0.7
  rho1 <- 0.5
  n=50
  K=7
 
  sa <- s(K, true_rho, true_rho1, n)
  
 
  mle_result <- optim(par = c(rho, rho1), fn = ll, hessian = TRUE, K = K, n = n, sa=sa)
  
 
  mle_rho <- mle_result$par[1]
  mle_rho1 <- mle_result$par[2]
  
  
  
  mse_rho <- MSE(mle_rho, true_rho)
  sd_rho <- sqrt(diag(solve(mle_result$hessian)))[1]
  
  cil_rho <- mle_rho - (1.96 * (sd_rho / sqrt(n)))
  ciu_rho <- mle_rho + (1.96 * (sd_rho / sqrt(n)))
  length_rho <- ciu_rho - cil_rho
  
  mle <- c(mle_rho,mse_rho,cil_rho,ciu_rho,length_rho)
  
  
  
  
  cat("True rho:", true_rho, "\n")
  cat("MLE of rho:", mle_rho, "\n")
  cat("MSE of rho:", mse_rho, "\n")
  
  mse_rho1 <- MSE(mle_rho1, true_rho1)
  sd_rho1 <- sqrt(diag(solve(mle_result$hessian)))[2]
  
  
  cil_rho1 <- mle_rho1 - (1.96 * (sd_rho1 / sqrt(n)))
  ciu_rho1 <- mle_rho1 + (1.96 * (sd_rho1 / sqrt(n)))
  length_rho1 <- ciu_rho1 - cil_rho1
  
  
  mle1 <- c(mle_rho1,mse_rho1,cil_rho1,ciu_rho1,length_rho1)
  
  cat("MLE of rho:", mle, "\n")
  cat("MLE of rho1:", mle1, "\n")


#Program_Table2:

library(extraDistr)
library(MLmetrics)
i <- function(K, rho,rho1, n) {
  # To fix the output, you can use set.seed followed by any random value, for example: set.seed(123) or set.seed(134567). 
  p <- 200
  P0 <- 0
  P0 <- P0 + (((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))
  P0 <- 1/P0
  Pn <- rep(0, p)
  N <- rep(0, p)
  c <- 1
  for (i in 0:p) {
    if (i < K) {
      P <- (rho1^i) * P0
    } else {
      P <- ((rho1^(K - 1)) * (rho^(i - K + 1))) * P0
    }
    Pn[c] <- P
    N[c] <- i
    c <- c + 1
  }
  Pn <- Pn/sum(Pn)
  max <- 50
  init <- 0
  final <- max
  T <- 5000
  Priori_rho <- rbeta(T,1,3) # beta
  Priori_rho1 <- rbetapr(T,1,3,1) # beta second kind
  while (init < n) {
    Ta <- min(final, n) - init
    R <- sample(N, Ta, replace = TRUE, prob = Pn)
    V_rho <- rep(0, T)
    V_rho1 <- rep(0, T)
    for (j in 1:T) {
      rho_j <- Priori_rho[j]
      for (r in 1:T) {
        rho1_r <- Priori_rho1[r]
        P0 <- 0
        P0 <- P0 + (((1 - rho1_r^K) / (1 - rho1_r)) + ((rho_j * rho1_r^(K - 1)) / (1 - rho_j)))
        P0 <- 1 / P0
        E1 <- 1
        for (i in 1:Ta) {
          i1 <- R[i]
          if (i1 < K) {
            E1a <- (rho1_r^i1) * P0
            E1 <- E1 * E1a
          } else {
            E2a <- ((rho1_r^(K - 1)) * (rho_j^(i1 - K + 1))) * P0
            E1 <- E1 * E2a
          }
        }
        V_rho[j] <- E1
        V_rho1[r] <- E1
      }
    }
    V_rho <- V_rho / sum(V_rho)
    V_rho1 <- V_rho1 / sum(V_rho1)
    Posteriori_rho <- sample(Priori_rho, T, replace = TRUE, prob = V_rho)
    Posteriori_rho1 <- sample(Priori_rho1, T, replace = TRUE, prob = V_rho1)
    init <- min(final, n)
    final <- final + max
    Priori_rho <- Posteriori_rho
    Priori_rho1 <- Posteriori_rho1
  }
  
  #for rho
  self <- mean(Posteriori_rho)
  risk.self <- mean((Posteriori_rho - self)^2)
  self.est <- c(self,MSE(self,rho))
  cdl <- sort(Posteriori_rho)[5.0*length(Posteriori_rho)/100]
  cdu<- sort(Posteriori_rho)[95.0*length(Posteriori_rho)/100]
  length <- cdu-cdl
  
  print(paste("self.est = ", self.est))
  print(paste("cdl = ", cdl))
  print(paste("cdu = ", cdu))
  print(paste("length = ", length))
  print(paste("risk.self = ", risk.self))
  
  #for rho1
  self1 <- mean(Posteriori_rho1)
  risk.self1 <- mean((Posteriori_rho1 - self1)^2)
  self.est1 <- c(self1,MSE(self1,rho1))
  cdl1 <- sort(Posteriori_rho1)[5.0*length(Posteriori_rho1)/100]
  cdu1 <- sort(Posteriori_rho1)[95.0*length(Posteriori_rho1)/100]
  length1 <- cdu1-cdl1
  
  print(paste("self.est1 = ", self.est1))
  print(paste("cd11 = ", cdl1))
  print(paste("cdu1 = ", cdu1))
  print(paste("length1 = ", length1))
  print(paste("risk.self1 = ", risk.self1))
  
}

i(K=3,rho=0.7,rho1=0.5,n=50)

#Program_Table3:
library(extraDistr)
library(MLmetrics)
library(Rcpp)
library(RGeode)
library(MASS)
i <- function(K, rho,rho1, n) {
  # To fix the output, you can use set.seed followed by any random value, for example: set.seed(123) or set.seed(134567). 
  p <- 200
  P0 <- 0
  P0 <- P0 + (((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))
  P0 <- 1/P0
  Pn <- rep(0, p)
  N <- rep(0, p)
  c <- 1
  for (i in 0:p) {
    if (i < K) {
      P <- (rho1^i) * P0
    } else {
      P <- ((rho1^(K - 1)) * (rho^(i - K + 1))) * P0
    }
    Pn[c] <- P
    N[c] <- i
    c <- c + 1
  }
  Pn <- Pn/sum(Pn)
  max <- 50
  init <- 0
  final <- max
  T <- 5000
  range <- c(0,1)
  Priori_rho <- rgammatr(T,1,3,range) 
  Priori_rho1 <- rbetapr(T,1,3,1) 
  while (init < n) {
    Ta <- min(final, n) - init
    R <- sample(N, Ta, replace = TRUE, prob = Pn)
    V_rho <- rep(0, T)
    V_rho1 <- rep(0, T)
    for (j in 1:T) {
      rho_j <- Priori_rho[j]
      for (r in 1:T) {
        rho1_r <- Priori_rho1[r]
        P0 <- 0
        P0 <- P0 + (((1 - rho1_r^K) / (1 - rho1_r)) + ((rho_j * rho1_r^(K - 1)) / (1 - rho_j)))
        P0 <- 1 / P0
        E1 <- 1
        for (i in 1:Ta) {
          i1 <- R[i]
          if (i1 < K) {
            E1a <- (rho1_r^i1) * P0
            E1 <- E1 * E1a
          } else {
            E2a <- ((rho1_r^(K - 1)) * (rho_j^(i1 - K + 1))) * P0
            E1 <- E1 * E2a
          }
        }
        V_rho[j] <- E1
        V_rho1[r] <- E1
      }
    }
    V_rho <- V_rho / sum(V_rho)
    V_rho1 <- V_rho1 / sum(V_rho1)
    Posteriori_rho <- sample(Priori_rho, T, replace = TRUE, prob = V_rho)
    Posteriori_rho1 <- sample(Priori_rho1, T, replace = TRUE, prob = V_rho1)
    init <- min(final, n)
    final <- final + max
    Priori_rho <- Posteriori_rho
    Priori_rho1 <- Posteriori_rho1
  }
  
  #for rho
  self <- mean(Posteriori_rho)
  risk.self <- mean((Posteriori_rho - self)^2)
  self.est <- c(self,MSE(self,rho))
  cdl <- sort(Posteriori_rho)[5.0*length(Posteriori_rho)/100]
  cdu<- sort(Posteriori_rho)[95.0*length(Posteriori_rho)/100]
  length <- cdu-cdl
  
  print(paste("self.est = ", self.est))
  print(paste("cdl = ", cdl))
  print(paste("cdu = ", cdu))
  print(paste("length = ", length))
  print(paste("risk.self = ", risk.self))
  
  #for rho1
  self1 <- mean(Posteriori_rho1)
  risk.self1 <- mean((Posteriori_rho1 - self1)^2)
  self.est1 <- c(self1,MSE(self1,rho1))
  cdl1 <- sort(Posteriori_rho1)[5.0*length(Posteriori_rho1)/100]
  cdu1 <- sort(Posteriori_rho1)[95.0*length(Posteriori_rho1)/100]
  length1 <- cdu1-cdl1
  
  print(paste("self.est1 = ", self.est1))
  print(paste("cd11 = ", cdl1))
  print(paste("cdu1 = ", cdu1))
  print(paste("length1 = ", length1))
  print(paste("risk.self1 = ", risk.self1))
  
}

i(K=3,rho=0.7,rho1=0.5,n=50)

#Program_Table4:
library(extraDistr)
library(MLmetrics)
i <- function(K, rho,rho1, n) {
  # To fix the output, you can use set.seed followed by any random value, for example: set.seed(123) or set.seed(134567). 
  p <- 200
  P0 <- 0
  P0 <- P0 + (((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))
  P0 <- 1/P0
  Pn <- rep(0, p)
  N <- rep(0, p)
  c <- 1
  for (i in 0:p) {
    if (i < K) {
      P <- (rho1^i) * P0
    } else {
      P <- ((rho1^(K - 1)) * (rho^(i - K + 1))) * P0
    }
    Pn[c] <- P
    N[c] <- i
    c <- c + 1
  }
  Pn <- Pn/sum(Pn)
  max <- 50
  init <- 0
  final <- max
  T <- 5000
  Priori_rho <- runif(T,0,1) 
  Priori_rho1 <- rbetapr(T,1,3,1) 
  while (init < n) {
    Ta <- min(final, n) - init
    R <- sample(N, Ta, replace = TRUE, prob = Pn)
    V_rho <- rep(0, T)
    V_rho1 <- rep(0, T)
    for (j in 1:T) {
      rho_j <- Priori_rho[j]
      for (r in 1:T) {
        rho1_r <- Priori_rho1[r]
        P0 <- 0
        P0 <- P0 + (((1 - rho1_r^K) / (1 - rho1_r)) + ((rho_j * rho1_r^(K - 1)) / (1 - rho_j)))
        P0 <- 1 / P0
        E1 <- 1
        for (i in 1:Ta) {
          i1 <- R[i]
          if (i1 < K) {
            E1a <- (rho1_r^i1) * P0
            E1 <- E1 * E1a
          } else {
            E2a <- ((rho1_r^(K - 1)) * (rho_j^(i1 - K + 1))) * P0
            E1 <- E1 * E2a
          }
        }
        V_rho[j] <- E1
        V_rho1[r] <- E1
      }
    }
    V_rho <- V_rho / sum(V_rho)
    V_rho1 <- V_rho1 / sum(V_rho1)
    Posteriori_rho <- sample(Priori_rho, T, replace = TRUE, prob = V_rho)
    Posteriori_rho1 <- sample(Priori_rho1, T, replace = TRUE, prob = V_rho1)
    init <- min(final, n)
    final <- final + max
    Priori_rho <- Posteriori_rho
    Priori_rho1 <- Posteriori_rho1
  }
  
  #for rho
  self <- mean(Posteriori_rho)
  risk.self <- mean((Posteriori_rho - self)^2)
  self.est <- c(self,MSE(self,rho))
  cdl <- sort(Posteriori_rho)[5.0*length(Posteriori_rho)/100]
  cdu<- sort(Posteriori_rho)[95.0*length(Posteriori_rho)/100]
  length <- cdu-cdl
  
  print(paste("self.est = ", self.est))
  print(paste("cdl = ", cdl))
  print(paste("cdu = ", cdu))
  print(paste("length = ", length))
  print(paste("risk.self = ", risk.self))
  
  #for rho1
  self1 <- mean(Posteriori_rho1)
  risk.self1 <- mean((Posteriori_rho1 - self1)^2)
  self.est1 <- c(self1,MSE(self1,rho1))
  cdl1 <- sort(Posteriori_rho1)[5.0*length(Posteriori_rho1)/100]
  cdu1 <- sort(Posteriori_rho1)[95.0*length(Posteriori_rho1)/100]
  length1 <- cdu1-cdl1
  
  print(paste("self.est1 = ", self.est1))
  print(paste("cd11 = ", cdl1))
  print(paste("cdu1 = ", cdu1))
  print(paste("length1 = ", length1))
  print(paste("risk.self1 = ", risk.self1))
  
}

i(K=3,rho=0.7,rho1=0.5,n=50)

#Program_Table5:
library(extraDistr)
library(MLmetrics)
i <- function(K, rho,rho1, n) {
 # To fix the output, you can use set.seed followed by any random value, for example: set.seed(123) or set.seed(134567). 
  p <- 200
  P0 <- 0
  P0 <- P0 + (((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))
  P0 <- 1/P0
  Pn <- rep(0, p)
  N <- rep(0, p)
  c <- 1
  for (i in 0:p) {
    if (i < K) {
      P <- (rho1^i) * P0
    } else {
      P <- ((rho1^(K - 1)) * (rho^(i - K + 1))) * P0
    }
    Pn[c] <- P
    N[c] <- i
    c <- c + 1
  }
  Pn <- Pn/sum(Pn)
  max <- 50
  init <- 0
  final <- max
  T <- 5000
  Priori_rho <- rbeta(T,1,3) 
  Priori_rho1 <- rgamma(T,1,3) 
  while (init < n) {
    Ta <- min(final, n) - init
    R <- sample(N, Ta, replace = TRUE, prob = Pn)
    V_rho <- rep(0, T)
    V_rho1 <- rep(0, T)
    for (j in 1:T) {
      rho_j <- Priori_rho[j]
      for (r in 1:T) {
        rho1_r <- Priori_rho1[r]
        P0 <- 0
        P0 <- P0 + (((1 - rho1_r^K) / (1 - rho1_r)) + ((rho_j * rho1_r^(K - 1)) / (1 - rho_j)))
        P0 <- 1 / P0
        E1 <- 1
        for (i in 1:Ta) {
          i1 <- R[i]
          if (i1 < K) {
            E1a <- (rho1_r^i1) * P0
            E1 <- E1 * E1a
          } else {
            E2a <- ((rho1_r^(K - 1)) * (rho_j^(i1 - K + 1))) * P0
            E1 <- E1 * E2a
          }
        }
        V_rho[j] <- E1
        V_rho1[r] <- E1
      }
    }
    V_rho <- V_rho / sum(V_rho)
    V_rho1 <- V_rho1 / sum(V_rho1)
    Posteriori_rho <- sample(Priori_rho, T, replace = TRUE, prob = V_rho)
    Posteriori_rho1 <- sample(Priori_rho1, T, replace = TRUE, prob = V_rho1)
    init <- min(final, n)
    final <- final + max
    Priori_rho <- Posteriori_rho
    Priori_rho1 <- Posteriori_rho1
  }
  
  #for rho
  self <- mean(Posteriori_rho)
  risk.self <- mean((Posteriori_rho - self)^2)
  self.est <- c(self,MSE(self,rho))
  cdl <- sort(Posteriori_rho)[5.0*length(Posteriori_rho)/100]
  cdu<- sort(Posteriori_rho)[95.0*length(Posteriori_rho)/100]
  length <- cdu-cdl
  
  print(paste("self.est = ", self.est))
  print(paste("cdl = ", cdl))
  print(paste("cdu = ", cdu))
  print(paste("length = ", length))
  print(paste("risk.self = ", risk.self))
  
  #for rho1
  self1 <- mean(Posteriori_rho1)
  risk.self1 <- mean((Posteriori_rho1 - self1)^2)
  self.est1 <- c(self1,MSE(self1,rho1))
  cdl1 <- sort(Posteriori_rho1)[5.0*length(Posteriori_rho1)/100]
  cdu1 <- sort(Posteriori_rho1)[95.0*length(Posteriori_rho1)/100]
  length1 <- cdu1-cdl1
  
  print(paste("self.est1 = ", self.est1))
  print(paste("cd11 = ", cdl1))
  print(paste("cdu1 = ", cdu1))
  print(paste("length1 = ", length1))
  print(paste("risk.self1 = ", risk.self1))
  
}

i(K=3,rho=0.7,rho1=0.5,n=50)


#Program_Table6:
library(extraDistr)
library(MLmetrics)
library(Rcpp)
library(RGeode)
library(MASS)
i <- function(K, rho,rho1, n) {
  # To fix the output, you can use set.seed followed by any random value, for example: set.seed(123) or set.seed(134567). 
  p <- 200
  P0 <- 0
  P0 <- P0 + (((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))
  P0 <- 1/P0
  Pn <- rep(0, p)
  N <- rep(0, p)
  c <- 1
  for (i in 0:p) {
    if (i < K) {
      P <- (rho1^i) * P0
    } else {
      P <- ((rho1^(K - 1)) * (rho^(i - K + 1))) * P0
    }
    Pn[c] <- P
    N[c] <- i
    c <- c + 1
  }
  Pn <- Pn/sum(Pn)
  max <- 50
  init <- 0
  final <- max
  T <- 5000
  range <- c(0,1)
  Priori_rho <- rgammatr(T,1,3,range) 
  Priori_rho1 <- rgamma(T,1,3) 
  while (init < n) {
    Ta <- min(final, n) - init
    R <- sample(N, Ta, replace = TRUE, prob = Pn)
    V_rho <- rep(0, T)
    V_rho1 <- rep(0, T)
    for (j in 1:T) {
      rho_j <- Priori_rho[j]
      for (r in 1:T) {
        rho1_r <- Priori_rho1[r]
        P0 <- 0
        P0 <- P0 + (((1 - rho1_r^K) / (1 - rho1_r)) + ((rho_j * rho1_r^(K - 1)) / (1 - rho_j)))
        P0 <- 1 / P0
        E1 <- 1
        for (i in 1:Ta) {
          i1 <- R[i]
          if (i1 < K) {
            E1a <- (rho1_r^i1) * P0
            E1 <- E1 * E1a
          } else {
            E2a <- ((rho1_r^(K - 1)) * (rho_j^(i1 - K + 1))) * P0
            E1 <- E1 * E2a
          }
        }
        V_rho[j] <- E1
        V_rho1[r] <- E1
      }
    }
    V_rho <- V_rho / sum(V_rho)
    V_rho1 <- V_rho1 / sum(V_rho1)
    Posteriori_rho <- sample(Priori_rho, T, replace = TRUE, prob = V_rho)
    Posteriori_rho1 <- sample(Priori_rho1, T, replace = TRUE, prob = V_rho1)
    init <- min(final, n)
    final <- final + max
    Priori_rho <- Posteriori_rho
    Priori_rho1 <- Posteriori_rho1
  }
  
  #for rho
  self <- mean(Posteriori_rho)
  risk.self <- mean((Posteriori_rho - self)^2)
  self.est <- c(self,MSE(self,rho))
  cdl <- sort(Posteriori_rho)[5.0*length(Posteriori_rho)/100]
  cdu<- sort(Posteriori_rho)[95.0*length(Posteriori_rho)/100]
  length <- cdu-cdl
  
  print(paste("self.est = ", self.est))
  print(paste("cdl = ", cdl))
  print(paste("cdu = ", cdu))
  print(paste("length = ", length))
  print(paste("risk.self = ", risk.self))
  
  #for rho1
  self1 <- mean(Posteriori_rho1)
  risk.self1 <- mean((Posteriori_rho1 - self1)^2)
  self.est1 <- c(self1,MSE(self1,rho1))
  cdl1 <- sort(Posteriori_rho1)[5.0*length(Posteriori_rho1)/100]
  cdu1 <- sort(Posteriori_rho1)[95.0*length(Posteriori_rho1)/100]
  length1 <- cdu1-cdl1
  
  print(paste("self.est1 = ", self.est1))
  print(paste("cd11 = ", cdl1))
  print(paste("cdu1 = ", cdu1))
  print(paste("length1 = ", length1))
  print(paste("risk.self1 = ", risk.self1))
  
}

i(K=3,rho=0.7,rho1=0.5,n=50)


#Program_Table7:
library(extraDistr)
library(MLmetrics)
i <- function(K, rho,rho1, n) {
  # To fix the output, you can use set.seed followed by any random value, for example: set.seed(123) or set.seed(134567). 
  p <- 200
  P0 <- 0
  P0 <- P0 + (((1 - rho1^K) / (1 - rho1)) + ((rho * rho1^(K - 1)) / (1 - rho)))
  P0 <- 1/P0
  Pn <- rep(0, p)
  N <- rep(0, p)
  c <- 1
  for (i in 0:p) {
    if (i < K) {
      P <- (rho1^i) * P0
    } else {
      P <- ((rho1^(K - 1)) * (rho^(i - K + 1))) * P0
    }
    Pn[c] <- P
    N[c] <- i
    c <- c + 1
  }
  Pn <- Pn/sum(Pn)
  max <- 50
  init <- 0
  final <- max
  T <- 5000
  Priori_rho <- runif(T,0,1) 
  Priori_rho1 <- rgamma(T,1,3) 
  while (init < n) {
    Ta <- min(final, n) - init
    R <- sample(N, Ta, replace = TRUE, prob = Pn)
    V_rho <- rep(0, T)
    V_rho1 <- rep(0, T)
    for (j in 1:T) {
      rho_j <- Priori_rho[j]
      for (r in 1:T) {
        rho1_r <- Priori_rho1[r]
        P0 <- 0
        P0 <- P0 + (((1 - rho1_r^K) / (1 - rho1_r)) + ((rho_j * rho1_r^(K - 1)) / (1 - rho_j)))
        P0 <- 1 / P0
        E1 <- 1
        for (i in 1:Ta) {
          i1 <- R[i]
          if (i1 < K) {
            E1a <- (rho1_r^i1) * P0
            E1 <- E1 * E1a
          } else {
            E2a <- ((rho1_r^(K - 1)) * (rho_j^(i1 - K + 1))) * P0
            E1 <- E1 * E2a
          }
        }
        V_rho[j] <- E1
        V_rho1[r] <- E1
      }
    }
    V_rho <- V_rho / sum(V_rho)
    V_rho1 <- V_rho1 / sum(V_rho1)
    Posteriori_rho <- sample(Priori_rho, T, replace = TRUE, prob = V_rho)
    Posteriori_rho1 <- sample(Priori_rho1, T, replace = TRUE, prob = V_rho1)
    init <- min(final, n)
    final <- final + max
    Priori_rho <- Posteriori_rho
    Priori_rho1 <- Posteriori_rho1
  }
  
  #for rho
  self <- mean(Posteriori_rho)
  risk.self <- mean((Posteriori_rho - self)^2)
  self.est <- c(self,MSE(self,rho))
  cdl <- sort(Posteriori_rho)[5.0*length(Posteriori_rho)/100]
  cdu<- sort(Posteriori_rho)[95.0*length(Posteriori_rho)/100]
  length <- cdu-cdl
  
  print(paste("self.est = ", self.est))
  print(paste("cdl = ", cdl))
  print(paste("cdu = ", cdu))
  print(paste("length = ", length))
  print(paste("risk.self = ", risk.self))
  
  #for rho1
  self1 <- mean(Posteriori_rho1)
  risk.self1 <- mean((Posteriori_rho1 - self1)^2)
  self.est1 <- c(self1,MSE(self1,rho1))
  cdl1 <- sort(Posteriori_rho1)[5.0*length(Posteriori_rho1)/100]
  cdu1 <- sort(Posteriori_rho1)[95.0*length(Posteriori_rho1)/100]
  length1 <- cdu1-cdl1
  
  print(paste("self.est1 = ", self.est1))
  print(paste("cd11 = ", cdl1))
  print(paste("cdu1 = ", cdu1))
  print(paste("length1 = ", length1))
  print(paste("risk.self1 = ", risk.self1))
  
}

i(K=3,rho=0.7,rho1=0.5,n=50)
