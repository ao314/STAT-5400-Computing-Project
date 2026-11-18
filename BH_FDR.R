FDR <- function(p_vals, truths, alpha=.05) {
  id <- order(p_vals)
  pvals <- sort(p_vals)
  truths <- as.logical(truths[id]) 
  m <- length(pvals)
  th <- (seq_len(m) / m) * alpha
  passed <- which(pvals <= th)
  
  if (length(passed) == 0) {
    k <- 0
  }
  else {
    k <- max(passed)
  }
  
  rejected <- seq_len(m) <= k
  
  null_rej <- sum(truths & rejected) # FP
  null_act <- sum(truths & !rejected) # TN
  alt_rej <- sum(!truths & rejected) # TP
  alt_act <- sum(!truths & !rejected) # FN
  
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