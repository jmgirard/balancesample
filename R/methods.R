# Balancing methods. See registry.R for the contract.

# anticlust -------------------------------------------------------------------

check_anticlust_options <- function(
  moments = 2L,
  covariances = FALSE,
  k_neighbours = Inf,
  local_maximum = TRUE,
  repetitions = 1L
) {
  if (!is_whole_number(moments) || moments < 1 || moments > 4) {
    cli::cli_abort("{.arg moments} must be 1, 2, 3, or 4.", call = NULL)
  }
  if (!is_flag(covariances)) {
    cli::cli_abort(
      "{.arg covariances} must be {.code TRUE} or {.code FALSE}.",
      call = NULL
    )
  }
  if (
    !is_number(k_neighbours) ||
      k_neighbours < 1 ||
      (is.finite(k_neighbours) && !is_whole_number(k_neighbours))
  ) {
    cli::cli_abort(
      "{.arg k_neighbours} must be a whole number of at least 1, or {.code Inf}.",
      call = NULL
    )
  }
  if (!is_flag(local_maximum)) {
    cli::cli_abort(
      "{.arg local_maximum} must be {.code TRUE} or {.code FALSE}.",
      call = NULL
    )
  }
  if (!is_whole_number(repetitions) || repetitions < 1) {
    cli::cli_abort(
      "{.arg repetitions} must be a whole number of at least 1.",
      call = NULL
    )
  }
  invisible()
}

# k-plus anticlustering (Papenberg, 2024) with anticlust::fast_anticlustering():
# the k-means criterion on the balance columns plus their higher moments (and,
# optionally, pairwise products), so that set means, variances, and so on are
# similar. The exchange step only swaps rows that share a category, so a
# stratified random start keeps the exact counts from allocate_strata(). With
# `local_maximum = TRUE`, exchange passes repeat until a pass makes no swap.
# Each repetition starts from a new stratified assignment. Groups are handled
# by aggregation: anticlust's `must_link` needs fold sizes that groups can
# fill exactly, which real grouped data rarely allows.
method_anticlust <- function(
  x_scaled,
  sizes,
  exact = NULL,
  group = NULL,
  moments = 2L,
  covariances = FALSE,
  k_neighbours = Inf,
  local_maximum = TRUE,
  repetitions = 1L
) {
  check_anticlust_options(
    moments,
    covariances,
    k_neighbours,
    local_maximum,
    repetitions
  )

  k <- length(sizes)
  # With one row per set, an exchange can only relabel the sets.
  if (all(sizes == 1L)) {
    return(stratified_random(sizes, exact))
  }
  features <- kplus_features(x_scaled, moments, covariances)
  categories <- if (is.null(exact)) NULL else as.integer(exact)
  stratum_counts <- function(a) {
    table(exact, factor(a, levels = seq_len(k)), dnn = NULL)
  }

  # A data frame is never mistaken for a distance matrix.
  features_df <- as.data.frame(features)

  # An exchange pass only makes swaps that improve the objective, so a pass
  # with no swaps means a local maximum. The pass limit guards against
  # cycling from rounding error.
  run <- function() {
    init <- stratified_random(sizes, exact)
    a <- init
    for (pass in seq_len(if (local_maximum) 100L else 1L)) {
      b <- as.integer(anticlust::fast_anticlustering(
        features_df,
        K = a,
        k_neighbours = k_neighbours,
        categories = categories
      ))
      done <- identical(b, a)
      a <- b
      if (done) {
        break
      }
    }
    obj <- anticlust::variance_objective(features, a)
    if (
      !is.null(exact) && !identical(stratum_counts(a), stratum_counts(init))
    ) {
      cli::cli_abort(
        "Internal error: anticlust changed the exact counts.",
        .internal = TRUE
      )
    }
    list(assign = a, objective = obj)
  }

  best <- run()
  for (r in seq_len(repetitions - 1)) {
    res <- run()
    if (res$objective > best$objective) {
      best <- res
    }
  }
  best$assign
}

# Balance columns, their moments 2 to `moments` (squared, cubed, ... deviations
# from the mean), and, if `covariances`, the products of each pair of columns,
# all standardized. Constant columns are dropped. The powers are built here
# because anticlust::kplus_moment_variables() raises its growing output, not
# the input, to each power in anticlust 0.8.18.
kplus_features <- function(x, moments, covariances) {
  centered <- sweep(x, 2, colMeans(x), "-")
  f <- x
  for (p in seq_len(moments)[-1]) {
    f <- cbind(f, centered^p)
  }
  if (covariances && ncol(x) > 1) {
    pairs <- which(upper.tri(diag(ncol(x))), arr.ind = TRUE)
    f <- cbind(f, x[, pairs[, 1], drop = FALSE] * x[, pairs[, 2], drop = FALSE])
  }
  # The features are computed from standardized columns, so rounding error is
  # relative to their size, not to their spread.
  f <- unname(as.matrix(f))
  standardize_matrix(f, apply(abs(f), 2, max))
}

# twinning --------------------------------------------------------------------

# `parts` is set internally, not by users. twinning's native code crashes R
# for a `leaf_size` below 1 or NA, so it is checked here.
check_twinning_options <- function(
  strategy = 1L,
  leaf_size = 8L,
  parts = NULL
) {
  if (!is_whole_number(strategy) || !strategy %in% 1:3) {
    cli::cli_abort("{.arg strategy} must be 1, 2, or 3.", call = NULL)
  }
  if (!is_whole_number(leaf_size) || leaf_size < 1) {
    cli::cli_abort(
      "{.arg leaf_size} must be a whole number of at least 1.",
      call = NULL
    )
  }
  invisible()
}

method_twinning <- function(
  x_scaled,
  sizes,
  exact = NULL,
  group = NULL,
  strategy = 1L,
  leaf_size = 8L,
  parts = NULL
) {
  rlang::check_installed("twinning", reason = "to use `method = \"twinning\"`.")
  check_twinning_options(strategy, leaf_size)
  k <- length(sizes)
  # twinning drops matrix dimensions for one column; a zero column does not
  # change any distance.
  if (ncol(x_scaled) == 1) {
    x_scaled <- cbind(x_scaled, 0)
  }

  # `parts` (set by balanced_initial_split() and
  # balanced_initial_validation_split()) gives the number of parts for each
  # set. Folds get one part each.
  blocks <- parts %||% rep.int(1L, k)
  m <- sum(blocks)

  # For two sets, strategies 1 and 2 of multiplet() reduce to twin(r = m),
  # which splits off one of `m` parts. Strategy 3 builds a different
  # sequence, so it uses multiplet().
  if (k == 2 && strategy != 3) {
    small <- which.min(sizes)
    other <- 3L - small
    split_fun <- function(x, sz) {
      idx <- twinning::twin(
        x,
        r = m,
        format_data = FALSE,
        leaf_size = leaf_size
      )
      a <- rep.int(other, nrow(x))
      a[idx] <- small
      fix_sizes(a, sz)
    }
    return(within_strata(x_scaled, sizes, exact, split_fun, min_rows = 2 * m))
  }

  if (strategy == 2 && !is_power_of_two(m)) {
    if (m == k) {
      cli::cli_abort(
        c(
          "Twinning strategy 2 needs {.arg v} to be a power of 2.",
          x = "{.arg v} is {k}.",
          i = "Use {.arg v} = 2, 4, 8, or 16, or use strategy 1 or 3."
        ),
        call = NULL
      )
    }
    cli::cli_abort(
      c(
        "Twinning strategy 2 needs each proportion to be a multiple of 1/m, where m is a power of 2.",
        x = "These proportions need m = {m}.",
        i = "Use proportions in quarters or eighths, such as {.code c(0.5, 0.25)}, or use strategy 1 or 3."
      ),
      call = NULL
    )
  }
  fold_fun <- function(x, sz) {
    a <- twinning::multiplet(
      x,
      k = m,
      strategy = strategy,
      format_data = FALSE,
      leaf_size = leaf_size
    )
    if (m == k) {
      return(relabel_to_sizes(a, sz))
    }
    target <- shuffle(rep.int(seq_len(k), blocks))
    fix_sizes(target[a], sz)
  }
  within_strata(x_scaled, sizes, exact, fold_fun, min_rows = 2 * m)
}

# Twinning parts for a two-way split with training proportion `prop`. r is
# the whole number nearest 1 / min(prop, 1 - prop), and `prop` must be within
# 1e-8 times 1 - 1/r (or 1/r, below 1/2) of that fraction. The comparison is
# with `prop` itself, because 1 - prop adds rounding error. The tolerance
# allows for proportions computed by subtraction, such as 1 - 0.9875.
# Returns c(r - 1, 1) or c(1, r - 1), or NULL. Like twin_parts(), it depends
# on `prop` only, not on the data size.
twin_split_parts <- function(prop) {
  r <- round(1 / min(prop, 1 - prop))
  if (!is.finite(r)) {
    return(NULL)
  }
  target <- if (prop >= 0.5) 1 - 1 / r else 1 / r
  if (abs(prop - target) > 1e-8 * target) {
    return(NULL)
  }
  if (prop >= 0.5) c(r - 1, 1) else c(1, r - 1)
}

# Number of twinning parts for each set: the smallest whole number `m` up to
# `max_parts` for which every proportion is a multiple of 1/m. Returns NULL if
# there is none. It depends on the proportions only, not on the data size.
twin_parts <- function(props, max_parts = 20L) {
  for (m in seq_len(max_parts)) {
    b <- props * m
    if (all(abs(b - round(b)) < 1e-8) && all(round(b) >= 1)) {
      return(as.integer(round(b)))
    }
  }
  NULL
}

is_power_of_two <- function(k) {
  k >= 1 && bitwAnd(as.integer(k), as.integer(k) - 1L) == 0L
}

# rerandomize -----------------------------------------------------------------

check_rerandomize_options <- function(p_accept = 0.1, max_draws = 1000L) {
  if (!is_number(p_accept) || p_accept <= 0 || p_accept > 1) {
    cli::cli_abort(
      "{.arg p_accept} must be a number above 0 and at most 1.",
      call = NULL
    )
  }
  if (!is_whole_number(max_draws) || max_draws < 1) {
    cli::cli_abort(
      "{.arg max_draws} must be a whole number of at least 1.",
      call = NULL
    )
  }
  invisible()
}

# Morgan and Rubin (2012): draw stratified random splits and keep the first
# whose Mahalanobis imbalance is under a threshold. For more than two sets the
# statistic is the largest value over "one set versus the rest" comparisons.
# The threshold is the `p_accept` quantile of the statistic over pilot draws
# (at least 100, more for small `p_accept`, up to `max_draws`) made under the
# same constraints, so it accounts for exact strata, unequal
# set sizes, and the dependence between the comparisons.
method_rerandomize <- function(
  x_scaled,
  sizes,
  exact = NULL,
  group = NULL,
  p_accept = 0.1,
  max_draws = 1000L
) {
  check_rerandomize_options(p_accept, max_draws)
  k <- length(sizes)
  inv <- pseudo_inverse(stats::cov(x_scaled))
  if (inv$rank == 0) {
    return(stratified_random(sizes, exact))
  }
  imbalance <- function(a) max_mahalanobis(x_scaled, a, k, inv$inverse)

  n_pilot <- max(100L, min(ceiling(10 / p_accept), max_draws))
  pilot <- vapply(
    seq_len(n_pilot),
    function(i) imbalance(stratified_random(sizes, exact)),
    numeric(1)
  )
  threshold <- stats::quantile(pilot, p_accept, names = FALSE, type = 1)
  # With discrete data many draws can tie at the threshold. Draws at the
  # threshold are accepted with the probability that brings the estimated
  # acceptance rate to `p_accept`.
  tol <- 1e-9 * max(1, abs(threshold))
  p_below <- mean(pilot < threshold - tol)
  p_at <- mean(abs(pilot - threshold) <= tol)
  p_tie <- min(1, (p_accept - p_below) / p_at)

  best <- NULL
  best_val <- Inf
  for (d in seq_len(max_draws)) {
    a <- stratified_random(sizes, exact)
    val <- imbalance(a)
    if (val < best_val) {
      best <- a
      best_val <- val
    }
    if (
      val < threshold - tol ||
        (val <= threshold + tol && stats::runif(1) < p_tie)
    ) {
      return(a)
    }
  }
  cli::cli_warn(c(
    "No draw met the rerandomization threshold after {max_draws} draw{?s}.",
    i = "Using the best draw (imbalance {signif(best_val, 3)}, threshold {signif(threshold, 3)}).",
    i = "Raise {.arg p_accept} or {.arg max_draws} in {.arg options}."
  ))
  best
}

max_mahalanobis <- function(x, assign, k, s_inv) {
  n <- nrow(x)
  sums <- rowsum(x, assign, reorder = TRUE)
  total <- colSums(x)
  nk <- tabulate(assign, k)
  vals <- vapply(
    seq_len(k),
    function(j) {
      d <- sums[j, ] / nk[j] - (total - sums[j, ]) / (n - nk[j])
      nk[j] * (n - nk[j]) / n * drop(crossprod(d, s_inv %*% d))
    },
    numeric(1)
  )
  max(vals)
}

pseudo_inverse <- function(m) {
  s <- svd(m)
  tol <- max(dim(m)) * max(s$d) * .Machine$double.eps
  keep <- s$d > tol
  if (!any(keep)) {
    return(list(inverse = m * 0, rank = 0L))
  }
  inverse <- s$v[, keep, drop = FALSE] %*%
    (t(s$u[, keep, drop = FALSE]) / s$d[keep])
  list(inverse = inverse, rank = sum(keep))
}

# interaction -----------------------------------------------------------------

check_interaction_options <- function(n_bins = 4L, pool = 0.1) {
  if (!is_whole_number(n_bins) || n_bins < 2) {
    cli::cli_abort(
      "{.arg n_bins} must be a whole number of at least 2.",
      call = NULL
    )
  }
  if (!is_number(pool) || pool < 0 || pool >= 1) {
    cli::cli_abort(
      "{.arg pool} must be a number from 0 up to (not including) 1.",
      call = NULL
    )
  }
  invisible()
}

# Baseline: cut each numeric column into quantile bins, cross all bins, pool
# sparse cells, and stratify at random. Exact strata stay exact because the
# crossed cells are nested inside them.
method_interaction <- function(
  x_scaled,
  sizes,
  exact = NULL,
  group = NULL,
  n_bins = 4L,
  pool = 0.1
) {
  check_interaction_options(n_bins, pool)
  bins <- lapply(seq_len(ncol(x_scaled)), function(j) {
    bin_column(x_scaled[, j], n_bins)
  })
  fine <- combine_levels(bins)
  fine <- pool_strata(fine, exact, min_size = pool * nrow(x_scaled))
  stratified_random(sizes, exact, fine)
}

# Quantile bins from ranks. Tied values share a rank, so they stay in one
# bin, and a heavily tied value no longer collapses all the bins into one. It
# can still share its bin with other values.
bin_column <- function(v, n_bins) {
  if (length(unique(v)) <= n_bins) {
    return(factor(v))
  }
  r <- rank(v, ties.method = "average")
  factor(ceiling(r * n_bins / length(v)))
}

pool_strata <- function(fine, exact, min_size) {
  prefix <- if (is.null(exact)) "" else paste0(as.integer(exact), "\r")
  key <- paste0(prefix, as.character(fine))
  tab <- table(key)
  small <- key %in% names(tab)[tab < min_size]
  key[small] <- paste0(prefix[pmin(length(prefix), which(small))], ".pooled")
  factor(key)
}

# option helpers ---------------------------------------------------------------

# One plain number that is not NA. Classed numbers (such as bit64's
# integer64) and arrays are refused, because their arithmetic can pass these
# checks while the native code of the methods reads something else.
# is.object() comes first, so that no class method of is.numeric() runs.
is_number <- function(x) {
  !is.object(x) &&
    !isS4(x) &&
    is.numeric(x) &&
    length(x) == 1 &&
    is.null(dim(x)) &&
    !is.na(x)
}

# TRUE or FALSE, as a plain vector of length one.
is_flag <- function(x) {
  !is.object(x) && !isS4(x) && rlang::is_bool(x) && is.null(dim(x))
}

# A whole number small enough to be an integer, so that seq_len() and the
# native code of the methods can use it.
is_whole_number <- function(x) {
  is_number(x) &&
    is.finite(x) &&
    x == round(x) &&
    abs(x) <= .Machine$integer.max
}

check_flag <- function(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!is_flag(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be {.code TRUE} or {.code FALSE}.",
      call = call
    )
  }
  invisible(x)
}
