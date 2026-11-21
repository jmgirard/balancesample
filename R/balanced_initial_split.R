#' Balanced training/testing split
#'
#' `balanced_initial_split()` makes one training/testing split in which the
#' two sets are similar on several variables at once. Use [rsample::training()]
#' and [rsample::testing()] to get the two data sets.
#'
#' @inheritParams balanced_vfold_cv
#' @param prop The proportion of data (or of groups, if you use `group`) for
#'   the training set. With `method = "twinning"`, `prop` must be `1 - 1/r`
#'   or `1/r` for a whole number `r`, for example 1/2, 2/3, 3/4, 9/10, or
#'   1/4. Write it as a fraction such as `2 / 3`, not a rounded decimal such
#'   as `0.67`. A value within a relative 1e-8 of such a fraction counts as
#'   that fraction.
#'
#' @details
#' The methods, options, and group handling are the same as for
#' [balanced_vfold_cv()]. With `method = "twinning"`, the split uses
#' [twinning::twin()], and the options are `leaf_size` only. Strata with
#' fewer than `2 * r` rows (or groups) are split at random.
#'
#' If a workflow has a postprocessor that must be trained, [tune::last_fit()]
#' splits the training set again with `balanced_initial_split()`, using the
#' same arguments.
#'
#' @return An `rsplit` object with class `balanced_initial_split` and
#'   `initial_split`.
#'
#' @seealso [balanced_vfold_cv()], [balance_check()]
#'
#' @examples
#' set.seed(1)
#' split <- balanced_initial_split(
#'   mtcars,
#'   prop = 3 / 4,
#'   balance = c(mpg, hp, wt),
#'   exact = am
#' )
#' split
#' rsample::training(split)
#' @export
balanced_initial_split <- function(
  data,
  prop = 3 / 4,
  balance = NULL,
  exact = NULL,
  group = NULL,
  method = "anticlust",
  options = list(),
  ...
) {
  call <- rlang::current_env()
  check_dots(..., call = call)
  if (
    !is.numeric(prop) ||
      length(prop) != 1 ||
      is.na(prop) ||
      prop <= 0 ||
      prop >= 1
  ) {
    cli::cli_abort("{.arg prop} must be a number between 0 and 1.", call = call)
  }
  entry <- get_balance_method(method, call = call)
  if (method == "twinning" && "strategy" %in% names(options)) {
    cli::cli_abort(
      "{.fn balanced_initial_split} doesn't use {.arg strategy}. Remove it from {.arg options}.",
      call = call
    )
  }
  check_method_options(entry, options, call = call)
  prep <- prepare_balance_data(
    data,
    rlang::enquo(balance),
    rlang::enquo(exact),
    rlang::enquo(group),
    call = call
  )
  if (!is.null(prep$group) && nlevels(prep$group) < 2) {
    cli::cli_abort("{.arg group} must have at least two groups.", call = call)
  }

  props <- c(prop, 1 - prop)
  method_options <- options
  if (method == "twinning") {
    method_options$parts <- twin_split_parts(prop)
    if (is.null(method_options$parts)) {
      cli::cli_abort(
        c(
          "{.code method = \"twinning\"} can only split off 1/r of the data, where r is a whole number.",
          i = "Use {.arg prop} = 1/2, 2/3, 3/4, 9/10, 1/4, or another value of 1 - 1/r or 1/r."
        ),
        call = call
      )
    }
    # A near miss such as 0.500000003 gives the same set sizes as 1/2.
    props <- method_options$parts / sum(method_options$parts)
  }
  a <- assign_balanced(prep, props, entry, method_options, call = call)
  warn_sparse_strata(prep, a, length(props), call = call)
  res <- rsample::make_splits(
    list(analysis = which(a == 1L), assessment = which(a == 2L)),
    data = data,
    class = c("balanced_initial_split", "initial_split")
  )

  attrib <- list(
    prop = prop,
    balance = prep$balance_names,
    exact = prep$exact_names,
    group = prep$group_name,
    method = method,
    options = options
  )
  for (nm in names(attrib)) {
    attr(res, nm) <- attrib[[nm]]
  }
  res
}

#' @export
print.balanced_initial_split <- function(x, ...) {
  details <- attributes(x)
  # The empty fallback of an inner calibration split stores no arguments.
  if (!is.null(details$method)) {
    cat_header(c(
      paste0("Balanced training/testing split (", details$method, ")"),
      design_parts(details)
    ))
  }
  NextMethod()
  invisible(x)
}
