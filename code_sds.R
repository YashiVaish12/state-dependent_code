#function for generate estimates
#SDD <- function(A, B, C, D) {
 # alapharaset1 <- 200
#X <-  (((1 - C^A) / (1 - C)) + ((B * C^(A - 1)) / (1 - B)))^-1
  #Xn <- rep(0, alapharaset1) N <- rep(0, alapharaset1)cc <- 1
  #for (i in 0:alapharaset1) {if (i < A) {P <- (C^i) * X} else {P <- ((C^(A - 1)) * (B^(i - A + 1))) * X}
   # Xn[cc] <- PN[cc]  icc <- cc + 1}
  #Xn <- Xn/sum(Xn);ms <- 50stt <- 0end <- ms
  #alapharaset <- 500 Alpha_B <- runif(alapharaset,0,1) Alpha_C <- rgamma(alapharaset,1,3) 
  #while (stt < D) {Ta <- min(end, D) - stt
#random sample generation
   # R <- sample(N, Ta, replace = TRUE, prob = Xn) Ve_B <- rep(0, alapharaset) Ve_C <- rep(0, alapharaset)
    #for (j in 1:alapharaset) {B_j <- Alpha_B[j]for (r in 1:alapharaset) {C_r <- Alpha_C[r] 
    #X <-  (((1 - C_r^A) / (1 - C_r)) + ((B_j * C_r^(A - 1)) / (1 - B_j)))^-1; E1 <- 1
        #for (i in 1:Ta) {i1 <- R[i] if (i1 < A) {E1a <- (C_r^i1) * X
            #E1 <- E1 * E1a} else {E2a <- ((C_r^(A - 1)) * (B_j^(i1 - A + 1))) * XE1 <- E1 * E2a}}
        #Ve_B[j] <- E1 Ve_C[r] <- E1}} Ve_B <- Ve_B / sum(Ve_B) Ve_C <- Ve_C / sum(Ve_C)
    #Pt_B <- sample(Alpha_B, alapharaset, replace = TRUE, prob = Ve_B )Pt_C <- sample(Alpha_C, alapharaset, replace = TRUE, prob = Ve_C)
    #stt <- min(end, D)end <- end + ms Alpha_B <- Pt_B Alpha_C <- Pt_C}
#bayes estimation
  #sf <- mean(Pt_B)print(paste("set = ", sf))
  #sf1 <- mean(Pt_C)print(paste("set1 = ", sf1))}

#SDD(A=3, B=0.7, C=0.5, D=50)
