
set.seed(5400)
total_tests <- c(4, 8, 16, 32, 64)
L <- c(5, 10)
ways <- c('d', 'c', 'i')
col <- c('a', 'b', 'c')
pnull <- c(0, 0.25, 0.5, 0.75)


for (let in c('a','b','c') )
{
    for(pnull in c(0, 0.25, 0.5, 0.75))
    {
        lname <-paste0('c', let, '_p', pnull )
        load(paste0("Data/", lname, ".RData"))
      }
}

for(rho in c(0.1, 0.3, 0.5)) {
    lname <- paste0('Data/dept_data_rho_', rho, '.RData')
    load(lname)
}


for(func in c("BY", "hochberg", "bonferroni", "FDR"))
{
    source(paste0(func, ".R"))
    for (let in c('a','b','c') )
        {
        for(pnull in c(0, 0.25, 0.5, 0.75))
            {
            lname <-paste0('c', let, '_p', pnull )
            assign(paste0(func, "_c", let, "_p", pnull),lapply(get(lname), function(pair) {
            A <- pair$p_vals
            B <- pair$truths
            mean(vapply(seq_len(nrow(A)),
                        function(i) do.call(paste0(func, "P"), list(A[i, ], B[i, ])), numeric(1)), na.rm = TRUE)
                }))
            }
        }
}

for(func in c("BY", "hochberg", "bonferroni", "FDR"))
{
    for(let in c('a','b','c'))
    {
        for(pnull in c(0, 0.25, 0.5, 0.75))
        {
            lname <- paste0(func, "_c", let, "_p", pnull)
            x <- unlist(get(lname))
            odd  <- seq(1, length(x), by = 2)    # positions 1, 3, 5, ...
            avgs <- (x[odd] + x[odd + 1]) / 2
            names(avgs) <- c("16", "32", "4", "64", "8")         # name each average after the first of its pair
            assign(lname, avgs)
        }
    }
}

for(let in c('a', 'b', 'c'))
{
    for(pnull in c(0, 0.25, 0.5, 0.75)) 
    {
        lname <- paste0("graph_c", let, "_p", pnull)
        assign(lname, list(
            BY = get(paste0("BY_c", let, "_p", pnull)),
            hochberg = get(paste0("hochberg_c", let, "_p", pnull)),
            bonferroni = get(paste0("bonferroni_c", let, "_p", pnull)),
            FDR = get(paste0("FDR_c", let, "_p", pnull))
        ))
        save(list = lname, file = paste0("Data/Graphs/", lname, ".RData"))
    }
}


# Dept Data

methods <- c("BY", "hochberg", "bonferroni", "FDR")   # must match your .R file names and <method>P functions

# source each method file once, not inside every iteration
for (meth in methods) source(paste0(meth, ".R"))

for (rho in c(0.1, 0.3, 0.5)) {

  # load the file into its own environment and grab the object inside,
  # so it doesn't matter what that object was called when it was saved
  e <- new.env()
  load(paste0("Data/dept_data_rho_", rho, ".RData"), envir = e)
  dat <- get(ls(e)[1], envir = e)       # list of lists, each with $p_vals and $truths

  res <- lapply(methods, function(meth) {
    fn <- get(paste0(meth, "P"))

    # average func over the rows of each pair of matrices
    avg_list <- lapply(dat, function(pair) {
      A <- pair$p_vals
      B <- pair$truths
      mean(vapply(seq_len(nrow(A)),
                  function(i) fn(A[i, ], B[i, ]),
                  numeric(1)), na.rm = TRUE)
    })

    # average every pair of consecutive numbers, named after the first of each pair
    x    <- unlist(avg_list)
    odd  <- seq(1, length(x), by = 2)
    avgs <- (x[odd] + x[odd + 1]) / 2
    names(avgs) <- c("16", "32", "4", "64", "8")
    avgs
  })
  names(res) <- methods

  lname <- paste0("graphs_dep_rho_", rho)
  assign(lname, res)
  save(list = lname, file = paste0("Data/Graphs/", lname, ".RData"))
}