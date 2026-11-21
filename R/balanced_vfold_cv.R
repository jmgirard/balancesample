#' Balanced V-fold cross-validation
#'
#' `balanced_vfold_cv()` splits the data into `v` folds of about equal size
#' that are similar on several variables at once. Each fold is used once as
#' the assessment set, and the other folds form the analysis set.
#'
#' @param data A data frame.
#' @param v The number of folds.
#' @param repeats The number of times to repeat the V-fold partitioning.
#' @param balance <[`tidy-select`][tidyselect::language]> Columns to balance
#'   across the sets. Numeric columns are standardized. So are `Date`,
#'   `POSIXct`, and `difftime` columns, as numbers of days or seconds. Factor,
#'   character, and logical columns become one indicator column per level.
#' @param exact <[`tidy-select`][tidyselect::language]> Categorical columns to
#'   stratify exactly. Each combination of levels is one stratum. In each
#'   stratum, the count in every set is the floor or the ceiling of its
#'   proportional share, based on the actual size of that set. For folds,
#'   the counts of a stratum also differ by at most one over the folds
#'   whenever both rules can hold at once. With `group`, these counts are
#'   numbers of groups, not rows. `exact` does the job of `strata` in
#'   rsample, but it takes several columns and its counts are exact. Numeric
#'   columns with whole-number values are treated as categories. Put
#'   continuous columns in `balance`, or use [cut()] to make a factor if you
#'   need exact counts within ranges. A stratum whose share of the smallest
#'   set is under one row (or group) can be missing from some sets. If more
#'   than half of the strata are missing from at least one set, you get a
#'   warning.
#' @param group <[`tidy-select`][tidyselect::language]> One column with group
#'   IDs. All rows of a group go to the same set. See Details.
#' @param method The balancing method: `"anticlust"`, `"twinning"`,
#'   `"rerandomize"`, or `"interaction"`. See Details.
#' @param options A named list of options for the method. See Details.
#' @param ... Not used. Must be empty.
#'
#' @details
#' # Methods
#'
#' * `"anticlust"` (default) uses k-plus anticlustering (Papenberg, 2024)
#'   with [anticlust::fast_anticlustering()]. It makes the means and the
#'   variances of the sets similar. Options are:
#'   * `moments` (default 2): the highest moment to balance. 1 balances
#'     means only, 3 adds skewness, and 4 adds kurtosis.
#'   * `covariances` (default `FALSE`): also balance the product of each
#'     pair of `balance` columns.
#'   * `k_neighbours` (default `Inf`): the number of nearest neighbours each
#'     row may swap with. `Inf` tries every swap, so the run time grows
#'     quickly with the number of rows. A small number, such as 20, is much
#'     faster on large data, but categorical `balance` columns can end up
#'     less balanced. Check the result with [balance_check()].
#'   * `local_maximum` (default `TRUE`): repeat the exchange passes until a
#'     pass makes no swap (at most 100 passes). `FALSE` runs one pass.
#'   * `repetitions` (default 1): the number of random starts. The best
#'     result is kept.
#' * `"twinning"` uses [twinning::multiplet()] for folds, or
#'   [twinning::twin()] for `v = 2` with strategy 1 or 2, which then give the
#'   same result. Options are `strategy` (1, 2, or 3; strategy 2 needs `v` to
#'   be a power of 2) and `leaf_size`. Strata with fewer than `2 * v` rows
#'   (or groups) are assigned at random. The twinning package is required.
#' * `"rerandomize"` draws stratified random splits and keeps the first one
#'   whose Mahalanobis imbalance is under a threshold (Morgan and Rubin,
#'   2012). For folds, the imbalance is the largest over the "one fold versus
#'   the rest" comparisons. The threshold is the `p_accept` quantile of the
#'   imbalance over pilot draws made with the same `exact` and `group`
#'   constraints (`max(100, min(ceiling(10 / p_accept), max_draws))`
#'   draws). Draws that tie with the threshold are accepted at random, so
#'   that the acceptance rate is close to `p_accept` even for discrete data.
#'   If `p_accept` is below `1 / max(100, max_draws)`, the pilot draws can't
#'   resolve it, and the threshold is the best pilot draw. Options are
#'   `p_accept`, the target acceptance rate (default 0.1), and `max_draws`
#'   (default 1000). If no draw meets the threshold, the best draw is used
#'   with a warning.
#' * `"interaction"` cuts numeric columns into quantile bins, crosses all
#'   bins, and stratifies at random. It is a baseline that shows how crossed
#'   strata become sparse. Options are `n_bins` (default 4) and `pool`
#'   (default 0.1): cells with less than this fraction of the rows are pooled.
#'
#' The method balances the columns left after standardizing: constant columns
#' are dropped, and with `group`, so are columns whose group means are all
#' equal up to rounding error (see Groups).
#' If none are left, for example with only `exact` columns and no `group`, the
#' method doesn't run and the sets are drawn at random, within the `exact`
#' strata if there are any. The limits the method itself checks, such as
#' twinning's strategy 2 needing a power of 2, then don't apply. The options
#' are still checked, and the split functions still check that twinning can
#' make `prop`.
#'
#' # Groups
#'
#' If you use `group`, each group becomes one row before the method runs.
#' This row holds the group means of the standardized `balance` columns and
#' the log of the group size. The group size is added only if a `balance`
#' column is left after standardizing or there are `exact` columns, and it is
#' dropped if all groups have the same size. The method assigns groups, and
#' then each row gets the set of its group. The number of groups per set
#' follows the target proportions, so the number of rows per set is only
#' approximate. `exact` columns must be constant within each group.
#'
#' # Use with tune
#'
#' tune and [rsample::reshuffle_rset()] find the function that made an
#' `rset` by its class name. They only find `balanced_vfold_cv()` when
#' balancesample is attached with `library(balancesample)`. If you call
#' `balancesample::balanced_vfold_cv()` without attaching the package,
#' [tune::fit_resamples()] and [tune::tune_grid()] fail with "Could not find
#' function".
#'
#' If a workflow has a postprocessor that must be trained, such as a
#' calibration from the tailor package, tune splits each analysis set again.
#' That inner split is a [balanced_initial_split()] with the same `balance`,
#' `exact`, `group`, `method`, and `options` (without twinning's
#' `strategy`), and with `prop = 1 - 1/v`.
#'
#' @return A tibble with class `balanced_vfold_cv` and `rset`. It has a
#'   `splits` column of `rsplit` objects and one or two ID columns.
#'
#' @references
#' Morgan, K. L., & Rubin, D. B. (2012). Rerandomization to improve covariate
#' balance in experiments. *The Annals of Statistics, 40*(2), 1263-1282.
#' \doi{10.1214/12-AOS1008}
#'
#' Papenberg, M. (2024). K-plus anticlustering: An improved k-means
#' criterion for maximizing between-group similarity. *British Journal of
#' Mathematical and Statistical Psychology, 77*(1), 80-102.
#' \doi{10.1111/bmsp.12315}
#'
#' Papenberg, M., & Klau, G. W. (2021). Using anticlustering to partition
#' data sets into equivalent parts. *Psychological Methods, 26*(2), 161-174.
#' \doi{10.1037/met0000301}
#'
#' Vakayil, A., & Joseph, V. R. (2022). Data twinning. *Statistical Analysis
#' and Data Mining, 15*(5), 598-610.
#' \doi{10.1002/sam.11574}
#'
#' @seealso [balanced_initial_split()], [balance_check()]
#'
#' @examples
#' set.seed(1)
#' folds <- balanced_vfold_cv(
#'   mtcars,
#'   v = 4,
#'   balance = c(mpg, hp, wt),
#'   exact = am
#' )
#' folds
#' balance_check(folds)
#'
#' # Other methods
#' balanced_vfold_cv(
#'   mtcars,
#'   v = 4,
#'   balance = c(mpg, hp),
#'   method = "rerandomize"
#' )
#' @examplesIf rlang::is_installed("twinning")
#' balanced_vfold_cv(mtcars, v = 4, balance = c(mpg, hp), method = "twinning")
#' @export
balanced_vfold_cv <- function(
  data,
  v = 10,
  repeats = 1,
  balance = NULL,
  exact = NULL,
  group = NULL,
  method = "anticlust",
  options = list(),
  ...
) {
  call <- rlang::current_env()
  check_dots(..., call = call)
  check_whole(v, min = 2, call = call)
  check_whole(repeats, min = 1, call = call)
  entry <- get_balance_method(method, call = call)
  check_method_options(entry, options, call = call)
  prep <- prepare_balance_data(
    data,
    rlang::enquo(balance),
    rlang::enquo(exact),
    rlang::enquo(group),
    call = call
  )

  n_units <- if (is.null(prep$group)) nrow(data) else nlevels(prep$group)
  unit <- if (is.null(prep$group)) "rows" else "groups"
  if (v > n_units) {
    cli::cli_abort(
      "{.arg v} ({v}) must not be more than the number of {unit} ({n_units}).",
      call = call
    )
  }

  props <- rep(1 / v, v)
  splits <- vector("list", v * repeats)
  for (r in seq_len(repeats)) {
    a <- assign_balanced(prep, props, entry, options, call = call)
    if (r == 1) {
      warn_sparse_strata(prep, a, v, call = call)
    }
    for (k in seq_len(v)) {
      splits[[(r - 1) * v + k]] <- rsample::make_splits(
        list(analysis = which(a != k), assessment = which(a == k)),
        data = data,
        class = "balanced_vfold_split"
      )
    }
  }

  if (repeats == 1) {
    ids <- tibble::tibble(id = names0(v, "Fold"))
  } else {
    ids <- tibble::tibble(
      id = rep(names0(repeats, "Repeat"), each = v),
      id2 = rep(names0(v, "Fold"), times = repeats)
    )
  }

  attrib <- list(
    v = v,
    repeats = repeats,
    balance = prep$balance_names,
    exact = prep$exact_names,
    group = prep$group_name,
    method = method,
    options = options
  )
  attrib <- attrib[
    !vapply(attrib, function(z) length(z) == 0 && !is.list(z), logical(1))
  ]

  rsample::new_rset(
    splits = splits,
    ids = ids,
    attrib = attrib,
    subclass = c("balanced_vfold_cv", "rset")
  )
}

#' @export
print.balanced_vfold_cv <- function(x, ...) {
  cat_header(header_parts(x))
  class(x) <- class(x)[!(class(x) %in% c("balanced_vfold_cv", "rset"))]
  print(x, ...)
  invisible(x)
}

#' @export
pretty.balanced_vfold_cv <- function(x, ...) {
  paste(header_parts(x, balance = FALSE), collapse = ", ")
}

# `balance = FALSE` leaves out the balance columns, so that pretty() stays
# short in the headers of other objects, such as nested_cv().
header_parts <- function(x, balance = TRUE) {
  details <- attributes(x)
  res <- paste0(
    details$v,
    "-fold balanced cross-validation (",
    details$method,
    ")"
  )
  if (!is.null(details$repeats) && details$repeats > 1) {
    res <- paste(res, "repeated", details$repeats, "times")
  }
  if (!balance) {
    details$balance <- NULL
  }
  c(res, design_parts(details))
}

# One part per column, so a long list can wrap between column names.
design_parts <- function(details) {
  res <- character()
  if (length(details$balance) > 0) {
    bal <- details$balance
    res <- c(res, paste("balancing", bal[1]), bal[-1])
  }
  if (length(details$exact) > 0) {
    exact <- details$exact
    res <- c(res, paste("stratified by", exact[1]), exact[-1])
  }
  if (length(details$group) > 0) {
    res <- c(res, paste("grouped by", details$group))
  }
  res
}

# Print a "# " header that joins `parts` with commas. A new line starts
# between parts when the next part, with its comma, would pass the console
# width. Widths are display widths, so wide characters count twice.
cat_header <- function(parts) {
  width <- getOption("width") - 2
  size <- function(x) nchar(x, type = "width")
  lines <- parts[1]
  n <- length(parts)
  for (i in seq_len(n)[-1]) {
    last <- length(lines)
    comma <- if (i < n) 1 else 0
    if (size(lines[last]) + 2 + size(parts[i]) + comma <= width) {
      lines[last] <- paste0(lines[last], ", ", parts[i])
    } else {
      lines[last] <- paste0(lines[last], ",")
      lines <- c(lines, parts[i])
    }
  }
  cat(paste0("# ", lines, "\n"), sep = "")
}

names0 <- function(num, prefix) {
  paste0(prefix, formatC(seq_len(num), width = nchar(num), flag = "0"))
}

check_whole <- function(
  x,
  min,
  max = Inf,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  ok <- is.numeric(x) &&
    length(x) == 1 &&
    is.finite(x) &&
    x == round(x) &&
    x >= min &&
    x <= max
  if (!ok) {
    cli::cli_abort(
      "{.arg {arg}} must be a whole number of at least {min}.",
      call = call
    )
  }
  invisible(x)
}

check_dots <- function(..., call = rlang::caller_env()) {
  if ("strata" %in% ...names()) {
    cli::cli_abort(
      c(
        "balancesample uses {.arg exact}, not {.arg strata}.",
        i = "Put categorical columns in {.arg exact}, for example {.code exact = class}.",
        i = "Put numeric columns in {.arg balance}."
      ),
      call = call
    )
  }
  # rsample arguments with a counterpart in balancesample.
  hints <- c(
    pool = "{.code method = \"interaction\"} pools small cells with {.code options = list(pool = 0.1)}.",
    breaks = "{.code method = \"interaction\"} sets the number of bins per numeric column with {.code options = list(n_bins = 4)}.",
    times = "{.fn balanced_vfold_cv} takes {.arg repeats}, which repeats the whole V-fold partition."
  )
  used <- intersect(names(hints), ...names())
  if (length(used) > 0) {
    cli::cli_abort(
      c(
        "balancesample doesn't use {.arg {used}}.",
        stats::setNames(hints[used], rep("i", length(used)))
      ),
      call = call
    )
  }
  rlang::check_dots_empty(call = call)
}
