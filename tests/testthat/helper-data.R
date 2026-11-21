sim_data <- function(n = 200, n_groups = 40) {
  withr::with_seed(42, {
    d <- data.frame(
      x1 = stats::rnorm(n),
      x2 = stats::rexp(n),
      x3 = stats::runif(n),
      cat = factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.6, 0.3, 0.1))),
      sex = sample(c("f", "m"), n, TRUE)
    )
    d$id <- sample(seq_len(n_groups), n, TRUE)
    # One site per group, so `site` can be stratified exactly with groups.
    d$site <- c("north", "south")[(d$id %% 2) + 1]
    d
  })
}

available_methods <- function() {
  m <- c("anticlust", "rerandomize", "interaction")
  if (rlang::is_installed("twinning")) {
    m <- c(m, "twinning")
  }
  m
}

# Each row is in exactly one assessment set per repeat.
assessed_once <- function(rset, n) {
  reps <- if ("id2" %in% names(rset)) rset$id else rep("r", nrow(rset))
  all(vapply(
    split(rset$splits, reps),
    function(sp) {
      out <- unlist(lapply(sp, rsample::complement))
      length(out) == n && setequal(out, seq_len(n))
    },
    logical(1)
  ))
}

# Counts per stratum and fold are within one of each other.
strata_exact <- function(rset, strata) {
  tab <- vapply(
    rset$splits,
    function(s) {
      tabulate(as.integer(strata)[rsample::complement(s)], nlevels(strata))
    },
    integer(nlevels(strata))
  )
  all(apply(tab, 1, function(r) max(r) - min(r) <= 1))
}

# Each stratum count in each assessment set is the floor or the ceiling of
# its share of that set's actual size.
shares_respected <- function(rset, strata) {
  m <- tabulate(as.integer(strata), nlevels(strata))
  n <- length(strata)
  all(vapply(
    rset$splits,
    function(s) {
      out <- rsample::complement(s)
      cnt <- tabulate(as.integer(strata)[out], nlevels(strata))
      share <- m * length(out) / n
      all(cnt >= floor(share + 1e-9) & cnt <= ceiling(share - 1e-9))
    },
    logical(1)
  ))
}

groups_intact <- function(rset, group) {
  all(vapply(
    rset$splits,
    function(s) {
      length(intersect(group[s$in_id], group[rsample::complement(s)])) == 0
    },
    logical(1)
  ))
}

max_abs_smd <- function(rset, ...) {
  bc <- balance_check(rset, ...)
  max(abs(bc$value[bc$statistic == "smd"]))
}
