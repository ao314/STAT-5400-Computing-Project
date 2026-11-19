    #Bonferroni controls Family-Wise Error Rate (probability of making any Type 1 errors in any comparisons)

#U <- True Null Hypothesis, Declared Non-Significant
#T <- Non-true Null Hypothesis, Declared Non-Significant
#V <- True Null Hypothesis, Declared Significant
#S <- Non-True Null Hypothesis, Declared Significant

#R <- V+S (observable random variable)
#m0 <- U+V
#m <- simulation repetitions

#PCER <- E(V/m)
#FWER <- P(V>=1)

#Q <- V/(V+S) (proportion of rejected null hypotheses that are erroneously rejected)
  #Q <- 0 if V+S=0

bonferroni <- function(pval, TRUEnull, a=0.05) {
  m  <- length(pval)               #hypotheses
  m0 <- sum(TRUEnull)              #true null hypotheses
  alpha <- a/m                     #adjusted alpha
  reject <- pval <= alpha          #rejected hypotheses
  V <- sum(TRUEnull & reject)      #true null , significant
  U <- sum(TRUEnull & !reject)     #true null, non-significant
  S <- sum(!TRUEnull & reject)     #non-true null, significant
  T <- sum(!TRUEnull & !reject)    #non-true null , non-significant
  power <- if(m-m0 <=0) NA_real_ else S/(m-m0)                #proportion of false null hypotheses, correctly rejected
  results <- data.frame(null_act=U, null_rej=V, alt_act=T, alt_rej=S, power=power)
  return(results)
}

#Just for testing:
#pval <- rnorm(16,mean=0.05,sd=0.1)
#TRUEnull <- c(1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,1)
#bon_method(pval,TRUEnull)
  




  
  
  
bonferroniP <- function(pval, TRUEnull, a=0.05) {
  m  <- length(pval)               #hypotheses
  m0 <- sum(TRUEnull)              #true null hypotheses
  alpha <- a/m                     #adjusted alpha
  reject <- pval <= alpha          #rejected hypotheses
  V <- sum(TRUEnull & reject)      #true null , significant
  U <- sum(TRUEnull & !reject)     #true null, non-significant
  S <- sum(!TRUEnull & reject)     #non-true null, significant
  T <- sum(!TRUEnull & !reject)    #non-true null , non-significant
  power <- S/(m-m0)                #proportion of false null hypotheses, correctly rejected
  return(power)
}

