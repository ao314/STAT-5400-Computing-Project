#' Joint Hypothesis Testing for p-values using the
#'  Hochberg (1988) procedure
#'
#' @param p_vals numeric vector containing all of the p-values
#' @param truths logical (or 0/1) vector: TRUE/1 if the null hypothesis is true,
#'   FALSE/0 if the alternative is true
#' @param alpha significance level (FWER level): default 0.05
#' 
#' @return
#' data.frame with 1 row and 5 columns:
#'      null_act (true nulls not rejected)
#'      null_rej (true nulls rejected, i.e. false discoveries)
#'      alt_act  (true alternatives not rejected, i.e. missed discoveries)
#'      alt_rej  (true alternatives rejected, i.e. true discoveries)
#'      power    (alt_rej / number of true alternatives; NA if there are none)



hochberg <- function(p_vals, truths, alpha= 0.05) {
  idx <- order(p_vals) # Get sorting order
  p_vals <- p_vals[idx] # sort p-values
  truths <- as.logical(truths[idx]) # sort hypotheses
  
  m <- length(p_vals)

  thresholds <- alpha / (m - seq_len(m) + 1) # Hochberg threshold: alpha / (m - k +1)
  passed <- which(p_vals<= thresholds)
  
  if (length(passed) == 0 ) {
    k <- 0
  } else {
    k <- max(passed)
  }
  
  rejected <- seq_len(m) <= k
  
  null_rej <- sum(truths & rejected)
  null_act <- sum(truths & !rejected)
  alt_rej <- sum(!truths & rejected)
  alt_act <- sum(!truths & !rejected)
  
  num_alt_hyp <- sum(!truths)
  
  if (num_alt_hyp != 0) {
    pwr <- alt_rej / num_alt_hyp
  } else {
    pwr <- NA_real_
  }
  
  
  out <- data.frame(null_act = null_act,
                    null_rej= null_rej, 
                    alt_act= alt_act,
                    alt_rej= alt_rej, 
                    power = pwr)
  return(out)
  
}