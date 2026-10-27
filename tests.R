# Method-agnostic test harness: compares your multiple-testing function
# against stats::p.adjust(). Your function must have the signature
#     fn(p_vals, truths, alpha)
# and return a 1-row data.frame with columns
#     null_act, null_rej, alt_act, alt_rej, power
#
# Usage:
#     source("BY.R"); source("test_methods.R")
#     test_method("BY")                         # looks up a function named BY
#     test_method("holm", fn = my_holm_func)    # or pass the function directly
#     test_method("BH", n_random = 2000, seed = 42)
#
# `method` must be one of p.adjust.methods:
#     "holm" "hochberg" "hommel" "bonferroni" "BH" "BY" "fdr" "none"

OUT_COLS <- c("null_act", "null_rej", "alt_act", "alt_rej", "power")

# Reference answer: reject H_i  <=>  adjusted p-value <= alpha
ref_output <- function(p, truths, alpha, method) {
  rej    <- p.adjust(p, method = method) <= alpha
  truths <- as.logical(truths)
  n_alt  <- sum(!truths)
  data.frame(null_act = sum(truths & !rej),
             null_rej = sum(truths & rej),
             alt_act  = sum(!truths & !rej),
             alt_rej  = sum(!truths & rej),
             power    = if (n_alt == 0) NA_real_ else sum(!truths & rej) / n_alt)
}

# Random dataset: nulls ~ Uniform(0,1), alternatives ~ Beta(0.1, 5)
random_data <- function(m, prop_null) {
  truths <- runif(m) < prop_null
  p <- ifelse(truths, runif(m), rbeta(m, 0.1, 5))
  list(p = p, truths = truths)
}

test_method <- function(method,
                        fn = get(method, mode = "function"),
                        n_random = 500,
                        seed = 1) {
  stopifnot(is.character(method), length(method) == 1,
            method %in% p.adjust.methods)
  set.seed(seed)

  results <- data.frame(test = character(), passed = logical(),
                        detail = character(), stringsAsFactors = FALSE)
  record <- function(test, passed, detail = "") {
    results[nrow(results) + 1, ] <<- list(test, isTRUE(passed), detail)
  }

  # Run fn, turning errors into failures
  run <- function(p, truths, alpha) {
    tryCatch(fn(p, truths, alpha), error = function(e) e)
  }

  # Compare one case against p.adjust; returns list(ok, detail)
  compare <- function(p, truths, alpha) {
    got <- run(p, truths, alpha)
    if (inherits(got, "error"))
      return(list(ok = FALSE, detail = paste("error:", conditionMessage(got))))
    exp <- ref_output(p, truths, alpha, method)
    ok  <- is.data.frame(got) && identical(names(got), OUT_COLS) &
           nrow(got) == 1 &&
           isTRUE(all.equal(got, exp, check.attributes = FALSE))
    detail <- if (ok) "" else
      paste0("got [", paste(unlist(got), collapse = ","), "] expected [",
             paste(unlist(exp), collapse = ","), "]")
    list(ok = ok, detail = detail)
  }

  # ---- 1. Output structure ----
  got <- run(c(.001, .5, .9), c(FALSE, TRUE, TRUE), .05)
  record("output is 1-row data.frame with correct columns",
         is.data.frame(got) && nrow(got) == 1 && identical(names(got), OUT_COLS))

  # ---- 2. Hand-picked edge cases (compared against p.adjust) ----
  edge <- list(
    "none rejected, all null"      = list(c(.9, .8, .7), rep(TRUE, 3), .05),
    "all rejected, all alternative"= list(rep(1e-9, 3), rep(FALSE, 3), .05),
    "all rejected, all null"       = list(rep(1e-9, 3), rep(TRUE, 3), .05),
    "none rejected, all alt"       = list(c(.9, .8, .7), rep(FALSE, 3), .05),
    "single hypothesis rejected"   = list(.03, FALSE, .05),
    "single hypothesis not rej."   = list(.07, FALSE, .05),
    "p-values of exactly 0 and 1"  = list(c(0, 1, 0, 1), c(FALSE, TRUE, FALSE, TRUE), .05),
    "tied p-values, mixed truths"  = list(rep(.01, 6), c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), .05),
    "unsorted input"               = list(c(.9, 1e-9, .8, 1e-9), c(TRUE, FALSE, TRUE, FALSE), .05),
    "first fails, later pass"      = list(c(.0095, .012, .02), rep(FALSE, 3), .05),
    "numeric 0/1 truths"           = list(c(1e-9, .9), c(0, 1), .05),
    "large m, all null"            = list(runif(1000), rep(TRUE, 1000), .05)
  )
  for (nm in names(edge)) {
    r <- compare(edge[[nm]][[1]], edge[[nm]][[2]], edge[[nm]][[3]])
    record(paste("edge:", nm), r$ok, r$detail)
  }

  # ---- 3. Random comparison against p.adjust ----
  fails <- 0; first_detail <- ""
  for (t in seq_len(n_random)) {
    m  <- sample(1:300, 1)
    d  <- random_data(m, prop_null = sample(c(0, .3, .7, .95, 1), 1))
    a  <- sample(c(.001, .01, .05, .1, .2), 1)
    r  <- compare(d$p, d$truths, a)
    if (!r$ok) { fails <- fails + 1; if (first_detail == "") first_detail <- r$detail }
  }
  record(sprintf("random: %d datasets match p.adjust", n_random),
         fails == 0,
         if (fails) sprintf("%d failures; first: %s", fails, first_detail) else "")

  # ---- 4. Property checks ----
  d <- random_data(100, .7); a <- .1
  o <- run(d$p, d$truths, a)
  record("counts sum to m",
         !inherits(o, "error") && sum(o[1:4]) == 100)

  perm <- sample(100)
  o2 <- run(d$p[perm], d$truths[perm], a)
  record("invariant to input order",
         !inherits(o, "error") && !inherits(o2, "error") &&
           isTRUE(all.equal(o, o2, check.attributes = FALSE)))

  nrej <- sapply(c(.01, .05, .1, .2, .5), function(al) {
    x <- run(d$p, d$truths, al)
    if (inherits(x, "error")) NA else x$null_rej + x$alt_rej
  })
  record("rejections nondecreasing in alpha",
         !anyNA(nrej) && all(diff(nrej) >= 0))

  pw <- run(d$p, d$truths, a)$power
  record("power in [0,1]", is.numeric(pw) && pw >= 0 && pw <= 1)

  # ---- Report ----
  cat(sprintf("\nTesting method '%s' against p.adjust\n", method))
  cat(strrep("-", 60), "\n")
  for (i in seq_len(nrow(results))) {
    cat(sprintf("%-4s %s\n", if (results$passed[i]) "PASS" else "FAIL",
                results$test[i]))
    if (!results$passed[i] && nzchar(results$detail[i]))
      cat("     ", results$detail[i], "\n")
  }
  cat(strrep("-", 60), "\n")
  cat(sprintf("%d / %d passed\n\n", sum(results$passed), nrow(results)))
  invisible(results)
}