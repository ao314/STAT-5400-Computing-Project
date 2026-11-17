Sim_Data <- function(m, p_null, div_way = c("c", "i", "d"), L, B = 20000, n = 1) {

    m_null <- p_null * m
    m_alt <- m - m_null
    truth <- c(rep(0, m_alt), rep(1, m_null))
    truths <- matrix(rep(truth, B), nrow= B, byrow = TRUE)
    

    mu_alt <- L * rep( c(1/4, 1/2, 3/4, 1), split_n(m_alt, div_way))
    mu_test <- c(mu_alt, rep(0, m_null))

    Z_stats <- sqrt(n) * matrix( rnorm( B * m, rep(mu_test, each = B), sd= 1/sqrt(n)), nrow= B)

    p_vals <- 2 * pnorm(- abs(Z_stats))

    return(list(p_vals = p_vals, truths = truths))
}




Sim_Data_Dep <- function(m, p_null, L, rho=0.1, B=20000) {
    Sigma <- matrix(rep(rho,m^2), nrow= m)
    diag(Sigma) <- rep(1, m)
    e <- eigen(Sigma)
    sqSigma <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
    m_null <- p_null * m
    m_alt <- m - m_null
    truth <- c(rep(0, m_alt), rep(1, m_null))
    truths <- matrix(rep(truth, B), nrow= B, byrow = TRUE)
    

    mu_alt <- L * rep( c(1/4, 1/2, 3/4, 1), split_n(m_alt, 'c'))
    mu_test <- c(mu_alt, rep(0, m_null))
    mu_mat <- matrix(mu_test, nrow = B, ncol = m, byrow = TRUE)
    Z_stats <- mu_mat + matrix( rnorm( B * m, 0, sd= 1), nrow= B) %*% sqSigma

    p_vals <- 2 * pnorm(- abs(Z_stats))

    return(list(p_vals = p_vals, truths = truths))
}



split_n <- function(n, type = c("c", "i", "d"), d = 1, k = 4) {
  type <- match.arg(type)
  stopifnot(n >= 0, k >= 1, d >= 0)

  idx <- 0:(k - 1)

  # constant split = step of 0
  if (type == "c") d <- 0

  # offsets from the base size: increasing 0..k-1, decreasing k-1..0
  offsets <- if (type == "d") rev(idx) else idx

  # if the step is too steep for n, reduce it until every group gets >= 1
  # (d = 0 is always allowed, so this always terminates)
  # d_req <- d
  while (d > 0 && (n - d * sum(idx)) %/% k < 1) d <- d - 1
  # if (d < d_req) warning("step reduced from ", d_req, " to ", d, " so n fits into ", k, " groups")

  # base size after removing the step pattern, rounded down
  a     <- (n - d * sum(idx)) %/% k
  sizes <- a + d * offsets

  # remainder (always < k): one extra unit per group
  r <- n - sum(sizes)
  if (r > 0) {
    pos <- if (type == "i") (k - r + 1):k else seq_len(r)
    sizes[pos] <- sizes[pos] + 1
  }

  sizes
}

set.seed(5400)
total_tests <- c(4, 8, 16, 32, 64)
L <- c(5, 10)
ways <- c('d', 'c', 'i')
col <- c('a', 'b', 'c')
pnull <- c(0, 0.25, 0.5, 0.75)


for ( m in total_tests ){
  for( l in L ) {
    for (w in 1:3) {
      for (p in pnull) {
        assign(paste0('m',m, '-l', l, '-c', col[w], '-p', p ), Sim_Data(m, p, ways[w], L = l))
      }
    }
  }
}
set.seed(5400)
for(rho in c(0.1, 0.3, 0.5)) {
  lname <- paste0('dept_data_rho_', rho)
  assign(lname, Sim_Data_Dep(m=16, p_null = 0.5, L = 5, rho = rho))
  save(list = lname , file = paste0("Data/", lname, ".RData"))
}



# save.image(file = "SimulationData.RData")
# load("SimulationData.RData")

# Find all objects matching the m-l-c-p pattern (separator can be _ or -)
objs <- ls(pattern = "^m[^_-]+[_-]l[^_-]+[_-]c[^_-]+[_-]p")

# Extract the "c-p" part: everything from the c-component onward
cp_ids <- sub("^m[^_-]+[_-]l[^_-]+[_-](c[^_-]+[_-]p.*)$", "\\1", objs)

for (cp in unique(cp_ids)) {
  # All objects sharing this c-p combination
  nms <- objs[cp_ids == cp]

  # Named list: names are the full object names (m-l-c-p)
  lst <- mget(nms)

  # Use a syntactically safe name for the saved object
  list_name <- gsub("-", "_", cp)
  assign(list_name, lst)
  save(list = list_name, file = paste0("Data/", list_name, ".RData"))
}
