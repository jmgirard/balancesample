# Methods for rsample generics that tune and other tidymodels packages call on
# resamples.

# tune calls internal_calibration_split() when a workflow has a postprocessor
# that must be trained (for example, a calibration from tailor). The analysis
# set is split again: the model is fit on one part and the postprocessor on the
# other. The inner split uses the same balancing as the outer one.

#' @exportS3Method rsample::internal_calibration_split
internal_calibration_split.balanced_vfold_split <- function(
  x,
  split_args,
  ...
) {
  rlang::check_dots_empty()
  if (is.null(split_args$prop)) {
    split_args$prop <- 1 - 1 / split_args$v
  }
  balanced_calibration_split(
    rsample::analysis(x),
    split_args,
    class = c(
      "balanced_vfold_split_cal",
      "internal_calibration_split",
      class(x)
    )
  )
}

#' @exportS3Method rsample::internal_calibration_split
internal_calibration_split.balanced_initial_split <- function(
  x,
  split_args,
  ...
) {
  rlang::check_dots_empty()
  # tune::last_fit() wraps the split in a manual_rset, so `split_args` is
  # empty there. The split stores its own arguments.
  stored <- attributes(x)[c(
    "prop",
    "balance",
    "exact",
    "group",
    "na",
    "method",
    "options"
  )]
  stored <- stored[!vapply(stored, is.null, logical(1))]
  split_args <- c(split_args, stored[setdiff(names(stored), names(split_args))])
  balanced_calibration_split(
    rsample::analysis(x),
    split_args,
    class = c(
      "balanced_initial_split_cal",
      "internal_calibration_split",
      class(x)
    )
  )
}

# As in rsample, the training set is split with the training share of the
# training and validation sets. twinning can only split off 1/r of the data,
# so with twinning the share is rounded to the nearest 1 - 1/r (or 1/r).
#' @exportS3Method rsample::internal_calibration_split
internal_calibration_split.balanced_initial_validation_split <- function(
  x,
  split_args,
  ...
) {
  rlang::check_dots_empty()
  stored <- attr(x, "val_att")
  stored <- stored[!vapply(stored, is.null, logical(1))]
  split_args <- c(split_args, stored[setdiff(names(stored), names(split_args))])
  balanced_calibration_split(
    rsample::training(x),
    split_args,
    class = c(
      "balanced_initial_validation_split_cal",
      "internal_calibration_split",
      "mc_split",
      "rsplit"
    ),
    share = three_way_fit_share
  )
}

# Errors here reach the fallback in balanced_calibration_split().
three_way_fit_share <- function(split_args) {
  p <- split_args$prop
  if (!is.numeric(p) || length(p) != 2 || !all(is.finite(p) & p > 0)) {
    cli::cli_abort(
      "{.arg prop} must be two finite numbers above 0.",
      call = NULL
    )
  }
  # Doubles, so that sums of large integers don't overflow.
  if (is.integer(p)) {
    p <- as.double(p)
  }
  twinning <- rlang::is_string(split_args$method, "twinning")
  if (twinning && sum(p) < 1) {
    # Use the whole-number parts of the outer split, so that near misses give
    # the same share as the fraction they match, and ties are exact.
    parts <- twin_parts(c(p, 1 - sum(p)))
    if (!is.null(parts)) {
      p <- parts[1:2]
    }
  }
  # Scale only if the sum overflows, so that other shares are unchanged.
  if (!is.finite(sum(p))) {
    p <- p / max(p)
  }
  if (twinning) {
    return(nearest_twin_share(p[[1]], p[[2]]))
  }
  p[[1]] / sum(p)
}

# The fit share nearest fit / (fit + cal) that twinning can make: 1 - 1/r, or
# 1/r if `fit` is the smaller. The shares 1/r are not evenly spaced, so the
# two nearest are compared. The comparison uses products, not quotients, so
# it is exact for the small whole-number parts of a twinning split. A tie
# goes to the larger r, which splits off the smaller set. A share that
# rounds to 1 is left for balanced_calibration_split() to reject.
nearest_twin_share <- function(fit, cal) {
  # Scale very large values so that the products can't overflow. Dividing by
  # a power of 2 is exact, so ties stay ties.
  if (max(fit, cal) > 1e150) {
    big <- 2^floor(log2(max(fit, cal)))
    fit <- fit / big
    cal <- cal / big
  }
  total <- fit + cal
  small <- min(fit, cal)
  r <- unique(pmax(2, c(floor(total / small), ceiling(total / small))))
  if (!all(is.finite(r))) {
    return(fit / total)
  }
  if (length(r) == 2) {
    # |1/r - small/total| is |total - r * small| / (r * total).
    near <- abs(total - r * small) * rev(r)
    r <- if (near[[2]] <= near[[1]]) r[[2]] else r[[1]]
  }
  if (fit >= cal) 1 - 1 / r else 1 / r
}

# Like rsample, fall back to an empty calibration set with a warning if the
# analysis set can't be split (for example, if it is too small).
# `share` turns `split_args` into the share of `analysis_set` to fit on.
balanced_calibration_split <- function(
  analysis_set,
  split_args,
  class,
  share = function(args) args$prop
) {
  # Warnings are held back and given once the split is made or has failed,
  # so that options(warn = 2) doesn't turn them into a failed split. They
  # come in their original order, before the fallback warning.
  warnings <- list()
  split <- tryCatch(
    withCallingHandlers(
      {
        options <- split_args$options %||% list()
        # The inner split is a balanced_initial_split(), which doesn't take
        # `strategy`. It also rejects options that aren't a list.
        if (is.list(options)) {
          options$strategy <- NULL
        }
        # A column can have no observed values in a part of the data even
        # when it has some in the whole. It carries no information there, so
        # it is left out rather than failing the split.
        balance <- drop_empty(analysis_set, split_args$balance)
        exact <- drop_empty(analysis_set, split_args$exact)
        balanced_initial_split(
          analysis_set,
          prop = share(split_args),
          balance = tidyselect::all_of(balance),
          exact = tidyselect::all_of(exact),
          group = tidyselect::all_of(split_args$group %||% character()),
          na = split_args$na %||% "error",
          method = split_args$method %||% "anticlust",
          options = options
        )
      },
      warning = function(cnd) {
        warnings[[length(warnings) + 1]] <<- cnd
        rlang::cnd_muffle(cnd)
      }
    ),
    error = identity
  )
  for (cnd in warnings) {
    warning(cnd)
  }
  if (inherits(split, "error")) {
    cli::cli_warn(
      "Cannot create calibration split; creating an empty calibration set.",
      parent = split
    )
    split <- rsample::make_splits(
      analysis_set,
      analysis_set[0, , drop = FALSE]
    )
  }
  class(split) <- class
  split
}

# The names in `cols` that `data` has, less those with no observed values.
drop_empty <- function(data, cols) {
  cols <- cols %||% character()
  empty <- vapply(
    cols,
    function(nm) nm %in% names(data) && all(is.na(data[[nm]])),
    logical(1)
  )
  cols[!empty]
}

#' @export
labels.balanced_vfold_cv <- function(object, make_factor = FALSE, ...) {
  if ("id2" %in% names(object)) {
    out <- paste(object$id, object$id2, sep = ".")
  } else {
    out <- as.character(object$id)
  }
  if (make_factor) {
    out <- as.factor(out)
  }
  out
}
