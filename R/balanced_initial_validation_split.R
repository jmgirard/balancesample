#' Balanced training/validation/testing split
#'
#' `balanced_initial_validation_split()` makes one three-way split in which
#' the training, validation, and testing sets are similar on several
#' variables at once. Use [rsample::training()], [rsample::validation()], and
#' [rsample::testing()] to get the three data sets, and
#' [rsample::validation_set()] to tune on the validation set.
#'
#' @inheritParams balanced_vfold_cv
#' @param prop A vector of two proportions: the first for the training set
#'   and the second for the validation set. The rest of the data goes to the
#'   testing set. With `group`, these are proportions of groups. With
#'   `method = "twinning"`, each of the three proportions must be a multiple
#'   of `1/m` for one whole number `m` up to 20, for example `c(0.6, 0.2)`
#'   (`m = 5`) or `c(0.5, 0.25)` (`m = 4`).
#'
#' @details
#' The three sets are made in one pass, so each one is balanced against the
#' other two. rsample's [rsample::initial_validation_split()] instead splits
#' off the training set first and then splits the rest.
#'
#' The methods, options, and group handling are the same as for
#' [balanced_vfold_cv()]. With `method = "twinning"`, the data are split
#' into `m` similar parts with [twinning::multiplet()], and each set gets a
#' whole number of parts. `strategy` 2 needs `m` to be a power of 2. Strata
#' with fewer than `2 * m` rows (or groups) are split at random.
#'
#' # Use with tune
#'
#' [rsample::validation_set()] turns the split into an `rset` for
#' [tune::tune_grid()] and [tune::fit_resamples()]. [tune::last_fit()] (tune
#' 1.1.2 or later) takes the split itself and fits on the training set, or on
#' the training and validation sets with `add_validation_set = TRUE`.
#'
#' If a workflow has a postprocessor that must be trained, tune splits the
#' training set again. Through [rsample::validation_set()], rsample makes
#' that inner split at random, without balancing. [tune::last_fit()] can't
#' train such a postprocessor on any three-way split, including those from
#' [rsample::initial_validation_split()] (as of tune 2.1.0).
#'
#' Calling [rsample::internal_calibration_split()] on this split yourself
#' gives a balanced split of the training set. As in rsample,
#' `prop[1] / sum(prop)` of the training set goes to the first part and the
#' rest to the calibration set. With `method = "twinning"`, this share is
#' rounded to the nearest `1 - 1/r` (or `1/r`, below 1/2) for a whole number
#' `r`, for example 14/17 to 5/6. As in rsample, if the training set can't
#' be split, you get a warning and an empty calibration set.
#'
#' @return A `three_way_split` object with class
#'   `balanced_initial_validation_split` and `initial_validation_split`.
#'
#' @seealso [balanced_initial_split()], [balance_check()]
#'
#' @examples
#' set.seed(1)
#' split <- balanced_initial_validation_split(
#'   mtcars,
#'   prop = c(0.6, 0.2),
#'   balance = c(mpg, hp, wt),
#'   exact = am
#' )
#' split
#' balance_check(split)
#' rsample::validation_set(split)
#' @export
balanced_initial_validation_split <- function(
  data,
  prop = c(0.6, 0.2),
  balance = NULL,
  exact = NULL,
  group = NULL,
  method = "anticlust",
  options = list(),
  na = "error",
  ...
) {
  call <- rlang::current_env()
  check_dots(..., call = call)
  if (
    !is.numeric(prop) ||
      length(prop) != 2 ||
      anyNA(prop) ||
      any(prop <= 0) ||
      sum(prop) >= 1
  ) {
    cli::cli_abort(
      c(
        "{.arg prop} must be two numbers above 0 that sum to less than 1.",
        i = "For example, {.code prop = c(0.6, 0.2)}."
      ),
      call = call
    )
  }
  entry <- get_balance_method(method, call = call)
  check_method_options(entry, options, call = call)
  prep <- prepare_balance_data(
    data,
    rlang::enquo(balance),
    rlang::enquo(exact),
    rlang::enquo(group),
    na = na,
    call = call
  )
  if (!is.null(prep$group) && nlevels(prep$group) < 3) {
    cli::cli_abort("{.arg group} must have at least three groups.", call = call)
  }

  props <- c(prop, 1 - sum(prop))
  method_options <- options
  if (method == "twinning") {
    method_options$parts <- twin_parts(props)
    if (is.null(method_options$parts)) {
      cli::cli_abort(
        c(
          "With {.code method = \"twinning\"}, each proportion must be a multiple of 1/m for one whole number m up to 20.",
          i = "For example, {.code prop = c(0.6, 0.2)} uses m = 5."
        ),
        call = call
      )
    }
    # Near misses give the same set sizes as the fractions they match.
    props <- method_options$parts / sum(method_options$parts)
  }
  a <- assign_balanced(prep, props, entry, method_options, call = call)
  warn_sparse_strata(prep, a, length(props), call = call)
  res <- list(
    data = data,
    train_id = which(a == 1L),
    val_id = which(a == 2L),
    test_id = NA,
    id = "split"
  )
  # rsample::validation_set() copies `val_att` to the attributes of the rset
  # it makes, so balance_check() finds the columns there too.
  attr(res, "val_att") <- list(
    prop = prop,
    balance = prep$balance_names,
    exact = prep$exact_names,
    group = prep$group_name,
    na = prep$na,
    method = method,
    options = options
  )
  class(res) <- c(
    "balanced_initial_validation_split",
    "initial_validation_split",
    "three_way_split"
  )
  res
}

#' @export
print.balanced_initial_validation_split <- function(x, ...) {
  details <- attr(x, "val_att")
  cat_header(c(
    paste0(
      "Balanced training/validation/testing split (",
      details$method,
      ")"
    ),
    design_parts(details)
  ))
  NextMethod()
  invisible(x)
}
