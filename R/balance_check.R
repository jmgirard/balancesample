#' Check the balance of a split or of a set of resamples
#'
#' `balance_check()` compares the assessment set of each split with its
#' analysis set. For a three-way split, it compares the validation set and the
#' testing set each with the training set. It reports a standardized mean
#' difference (SMD) for each variable and the energy distance for all
#' variables together. It works
#' with any `rset` or `rsplit`, so you can compare balanced folds with
#' ordinary [rsample::vfold_cv()] folds.
#'
#' @param x An `rset` (for example from [balanced_vfold_cv()]), an `rsplit`
#'   (for example from [balanced_initial_split()]), or a three-way split (from
#'   [balanced_initial_validation_split()] or
#'   [rsample::initial_validation_split()]).
#' @param ... For `balance_check()`, <[`tidy-select`][tidyselect::language]>
#'   columns to check. If empty, the `balance` and `exact` columns stored in
#'   `x` are used. For `summary()`, not used and must be empty.
#' @param max_n The largest number of rows per set for the energy distance.
#'   Larger sets are subsampled at random to keep memory use low.
#' @param variance_ratio If `TRUE`, also report the variance ratio of each
#'   numeric variable.
#'
#' @details
#' For a numeric variable, the SMD is the difference in means (assessment
#' minus analysis) divided by the pooled standard deviation,
#' `sqrt((var_assessment + var_analysis) / 2)`. For a categorical variable,
#' there is one SMD per level, from the proportions of that level. Values
#' under 0.1 in absolute size are often taken as good balance. If both sets
#' have zero variance but different means, the SMD is `Inf` or `-Inf`; the
#' plot shows these values at the edge of the panel.
#'
#' Stored `exact` columns with whole numbers are treated as categories, as
#' the split functions treat them, so they get one SMD per level and no
#' variance ratio. This holds whether you name them or not. Other
#' whole-number columns are numeric.
#'
#' With `variance_ratio = TRUE`, each numeric variable also gets the variance
#' of the assessment set divided by the variance of the analysis set. It is 1
#' when the spreads are equal, and values between 0.5 and 2 are often taken
#' as good balance. The SMD can be 0 while the spreads differ, and the
#' variance ratio shows which variable differs. A variance within rounding
#' error of zero counts as zero: if both are zero, the ratio is 1, and if one
#' is, it is 0 or `Inf`. With fewer than two values in a set, it is `NA`.
#'
#' Missing values are dropped variable by variable for the SMD and the
#' variance ratio. The energy
#' distance uses the rows that have no missing values in any checked column.
#'
#' The energy distance is the V-statistic
#' `2 E|X - Y| - E|X - X'| - E|Y - Y'|` on the standardized variables (one
#' indicator column per categorical level). It is 0 when the two sets have the
#' same distribution, and it also detects differences in spread and shape.
#'
#' Use `summary()` on the result to get one row per split.
#'
#' @return A tibble with class `balance_check`. It has the ID columns of `x`,
#'   and the columns `variable`, `level`, `statistic` (`"smd"`, `"energy"`, or
#'   `"variance_ratio"`), and `value`. Energy rows have `NA` in `variable` and
#'   `level`.
#'
#'   `summary()` returns a tibble with one row per split. It has the ID
#'   columns, `max_abs_smd` (the largest absolute SMD in the split), the
#'   `variable` and `level` with that SMD (the first one if several tie), and
#'   `energy`. With variance ratios, it also has `variance_ratio`, the ratio
#'   farthest from 1 on the log scale (so 0.25 is farther than 3). If the
#'   split has no SMD that isn't missing, `max_abs_smd`,
#'   `variable`, and `level` are `NA`. If it has no energy row, `energy` is
#'   `NA`. Each split must have its own ID, so no ID may have more than one
#'   energy row. Only the ID columns that `balance_check()` made are used, so
#'   columns you add are ignored. If you rename or drop an ID column, the
#'   columns whose names start with `id` are used instead.
#'
#' @examples
#' set.seed(1)
#' folds <- balanced_vfold_cv(mtcars, v = 4, balance = c(mpg, hp, wt))
#' balance_check(folds)
#' summary(balance_check(folds))
#' summary(balance_check(folds, variance_ratio = TRUE))
#'
#' # Compare with ordinary folds
#' summary(balance_check(rsample::vfold_cv(mtcars, v = 4), mpg, hp, wt))
#' @export
balance_check <- function(x, ..., max_n = 2000, variance_ratio = FALSE) {
  UseMethod("balance_check")
}

#' @export
balance_check.default <- function(
  x,
  ...,
  max_n = 2000,
  variance_ratio = FALSE
) {
  cli::cli_abort(
    "{.arg x} must be an {.cls rset}, an {.cls rsplit}, or a {.cls three_way_split}, not {.obj_type_friendly {x}}."
  )
}

#' @rdname balance_check
#' @export
balance_check.rset <- function(
  x,
  ...,
  max_n = 2000,
  variance_ratio = FALSE
) {
  check_whole(max_n, min = 2)
  check_flag(variance_ratio)
  if (nrow(x) == 0) {
    cli::cli_abort("{.arg x} has no splits.")
  }
  # The columns are selected in the first split. The splits of a
  # manual_rset() can have different data, so check that the others have them.
  cols <- check_columns(x$splits[[1]]$data, x, ...)
  for (i in seq_len(nrow(x))) {
    lacking <- setdiff(cols, names(x$splits[[i]]$data))
    if (length(lacking) > 0) {
      cli::cli_abort(
        c(
          "Every split of {.arg x} must have the columns to check.",
          x = "Split {i} has no {.field {lacking}}."
        )
      )
    }
  }
  id_cols <- grep("^id", names(x), value = TRUE)
  ids <- tibble::as_tibble(as.data.frame(x)[id_cols])

  res <- lapply(seq_len(nrow(x)), function(i) {
    # Splits of a manual_rset() can have different data.
    split <- x$splits[[i]]
    data <- exact_as_factors(split$data, x, cols)
    out <- split_balance(split, cols, max_n, data, variance_ratio)
    id_rows <- ids[rep(i, nrow(out)), , drop = FALSE]
    cbind(id_rows, out)
  })
  new_balance_check(do.call(rbind, res), id_cols)
}

#' @rdname balance_check
#' @export
balance_check.rsplit <- function(
  x,
  ...,
  max_n = 2000,
  variance_ratio = FALSE
) {
  check_whole(max_n, min = 2)
  check_flag(variance_ratio)
  cols <- check_columns(x$data, x, ...)
  out <- split_balance(
    x,
    cols,
    max_n,
    exact_as_factors(x$data, x, cols),
    variance_ratio
  )
  new_balance_check(cbind(tibble::tibble(id = rep("Split", nrow(out))), out))
}

#' @rdname balance_check
#' @export
balance_check.three_way_split <- function(
  x,
  ...,
  max_n = 2000,
  variance_ratio = FALSE
) {
  check_whole(max_n, min = 2)
  check_flag(variance_ratio)
  cols <- check_columns(x$data, x, ...)
  data <- exact_as_factors(x$data, x, cols)
  test_id <- setdiff(seq_len(nrow(x$data)), c(x$train_id, x$val_id))
  res <- lapply(
    list(Validation = x$val_id, Testing = test_id),
    function(out_id) {
      split <- rsample::make_splits(
        list(analysis = x$train_id, assessment = out_id),
        data = x$data
      )
      split_balance(split, cols, max_n, data, variance_ratio)
    }
  )
  ids <- rep(names(res), vapply(res, nrow, integer(1)))
  new_balance_check(cbind(tibble::tibble(id = ids), do.call(rbind, res)))
}

check_columns <- function(data, x, ..., call = rlang::caller_env()) {
  pos <- tidyselect::eval_select(
    rlang::expr(c(!!!rlang::enquos(...))),
    data,
    allow_rename = FALSE,
    error_call = call
  )
  cols <- names(pos)
  if (length(cols) == 0) {
    stored <- stored_args(x)
    cols <- unique(c(stored$balance, stored$exact))
  }
  if (length(cols) == 0) {
    cli::cli_abort(
      c(
        "Select the columns to check.",
        i = "For example, {.code balance_check(folds, age, income)}."
      ),
      call = call
    )
  }
  cols
}

# Three-way splits keep their arguments in `val_att`.
stored_args <- function(x) {
  attr(x, "val_att") %||% attributes(x)
}

# The split functions treat whole-number `exact` columns as categories, so
# the check does too, with one SMD per level.
exact_as_factors <- function(data, x, cols) {
  for (nm in intersect(cols, stored_args(x)$exact)) {
    if (is.numeric(data[[nm]])) {
      data[[nm]] <- factor(data[[nm]])
    }
  }
  data
}

split_balance <- function(
  split,
  cols,
  max_n,
  data = split$data,
  variance_ratio = FALSE
) {
  an <- split$in_id
  as <- rsample::complement(split)

  smd <- lapply(cols, function(nm) {
    v <- data[[nm]]
    if (is.logical(v)) {
      v <- factor(v, levels = c(FALSE, TRUE))
    }
    if (is.character(v)) {
      v <- factor(v)
    }
    if (is.factor(v)) {
      lv <- levels(droplevels(v))
      vals <- vapply(
        lv,
        function(l) {
          z <- as.numeric(v == l)
          smd_value(z[as], z[an], binary = TRUE)
        },
        numeric(1)
      )
      return(tibble::tibble(
        variable = nm,
        level = lv,
        statistic = "smd",
        value = unname(vals)
      ))
    }
    v <- as.numeric(as_balance_numeric(v))
    out <- tibble::tibble(
      variable = nm,
      level = NA_character_,
      statistic = "smd",
      value = smd_value(v[as], v[an])
    )
    if (variance_ratio) {
      out <- rbind(
        out,
        tibble::tibble(
          variable = nm,
          level = NA_character_,
          statistic = "variance_ratio",
          value = variance_ratio_value(v[as], v[an])
        )
      )
    }
    out
  })

  keep <- which(stats::complete.cases(data[cols]))
  pos <- integer(nrow(data))
  pos[keep] <- seq_along(keep)
  as_k <- pos[as][pos[as] > 0]
  an_k <- pos[an][pos[an] > 0]
  en <- NA_real_
  if (length(as_k) > 0 && length(an_k) > 0) {
    xs <- build_balance_matrix(data[keep, cols, drop = FALSE])
    en <- energy_distance(
      xs[as_k, , drop = FALSE],
      xs[an_k, , drop = FALSE],
      max_n
    )
  }
  energy <- tibble::tibble(
    variable = NA_character_,
    level = NA_character_,
    statistic = "energy",
    value = en
  )
  do.call(rbind, c(smd, list(energy)))
}

smd_value <- function(a, b, binary = FALSE) {
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]
  if (length(a) == 0 || length(b) == 0) {
    return(NA_real_)
  }
  if (binary) {
    pa <- mean(a)
    pb <- mean(b)
    diff <- pa - pb
    den <- sqrt((pa * (1 - pa) + pb * (1 - pb)) / 2)
    tol <- near_zero(1)
  } else {
    # Each variance is measured within its own set, before the shift below.
    den <- sqrt((set_variance(a) + set_variance(b)) / 2)
    # Measure from an observed value, so that a large origin (such as a
    # date-time in seconds) doesn't swamp small differences.
    ref <- c(a, b)[1]
    a <- a - ref
    b <- b - ref
    diff <- mean(a) - mean(b)
    tol <- near_zero(max(abs(c(a, b))))
  }
  if (den == 0) {
    return(if (abs(diff) <= tol) 0 else sign(diff) * Inf)
  }
  out <- diff / den
  if (abs(out) < 1e-12) 0 else out
}

# Variance of one set, or 0 if it is within rounding error of zero. Each set
# is measured from its own first value, so a gap between the sets doesn't
# hide the spread within them.
set_variance <- function(a) {
  if (length(a) < 2) {
    return(0)
  }
  a <- a - a[1]
  v <- stats::var(a)
  if (v <= near_zero(max(abs(a)))^2) 0 else v
}

# Variance of `a` over the variance of `b`. A variance within rounding error
# of zero counts as zero: two zeros give 1, and a zero over a variance above
# zero gives 0 (or Inf the other way round).
variance_ratio_value <- function(a, b) {
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]
  if (length(a) < 2 || length(b) < 2) {
    return(NA_real_)
  }
  va <- set_variance(a)
  vb <- set_variance(b)
  if (va == 0 && vb == 0) {
    return(1)
  }
  if (vb == 0) {
    return(Inf)
  }
  if (va == 0) {
    return(0)
  }
  va / vb
}

energy_distance <- function(a, b, max_n = 2000) {
  if (ncol(a) == 0) {
    return(0)
  }
  if (nrow(a) > max_n) {
    a <- a[sample.int(nrow(a), max_n), , drop = FALSE]
  }
  if (nrow(b) > max_n) {
    b <- b[sample.int(nrow(b), max_n), , drop = FALSE]
  }
  2 * mean_dist(a, b) - mean_dist(a, a) - mean_dist(b, b)
}

mean_dist <- function(a, b) {
  d2 <- outer(rowSums(a^2), rowSums(b^2), "+") - 2 * tcrossprod(a, b)
  mean(sqrt(pmax(d2, 0)))
}

# `id_cols` records the ID columns, so that summary() can tell them from
# columns added later.
new_balance_check <- function(x, id_cols = "id") {
  x <- tibble::as_tibble(x)
  attr(x, "id_cols") <- id_cols
  class(x) <- c("balance_check", class(x))
  x
}

#' @export
print.balance_check <- function(x, ...) {
  smd <- x$value[x$statistic == "smd"]
  en <- x$value[x$statistic == "energy"]
  n_split <- length(en)
  smd <- abs(smd[!is.na(smd)])
  max_smd <- if (length(smd) > 0) max(smd) else NA_real_
  mean_en <- if (any(!is.na(en))) mean(en, na.rm = TRUE) else NA_real_
  cat_header(c(
    paste0("Balance check: ", n_split, " split", if (n_split != 1) "s"),
    paste0("max |SMD| = ", format(signif(max_smd, 3))),
    paste0("mean energy distance = ", format(signif(mean_en, 3)))
  ))
  print(unclass_balance(x), ...)
  invisible(x)
}

#' @rdname balance_check
#' @param object A `balance_check` object.
#' @export
summary.balance_check <- function(object, ...) {
  rlang::check_dots_empty()
  x <- unclass_balance(object)
  stat_cols <- c("variable", "level", "statistic", "value")
  # Use the ID columns balance_check() recorded, so that added columns are
  # ignored. If the record is missing or an ID column was renamed or
  # dropped, take the columns whose names start with "id", as rsample does.
  id_cols <- attr(object, "id_cols")
  if (is.null(id_cols) || !all(id_cols %in% names(x))) {
    id_cols <- grep("^id", names(x), value = TRUE)
  }
  if (!all(stat_cols %in% names(x)) || length(id_cols) == 0) {
    cli::cli_abort(
      c(
        "{.arg object} must have ID columns and the columns {.code {stat_cols}}.",
        i = "Summarize the output of {.fn balance_check}."
      )
    )
  }
  stats <- c("smd", "energy", "variance_ratio")
  if (!all(x$statistic %in% stats)) {
    cli::cli_abort(
      "{.code statistic} must be {.or {.val {stats}}} in every row of {.arg object}."
    )
  }
  # Integer codes keep `NA` apart from "NA" and can't collide when pasted.
  codes <- lapply(x[id_cols], function(v) match(v, unique(v)))
  key <- do.call(paste, c(unname(codes), sep = "."))
  rows <- split(seq_len(nrow(x)), factor(key, levels = unique(key)))
  n_energy <- vapply(rows, function(i) sum(x$statistic[i] == "energy"), 1L)
  if (any(n_energy > 1)) {
    cli::cli_abort(
      c(
        "Each split in {.arg object} must have its own ID.",
        x = "Some IDs have more than one energy row."
      )
    )
  }
  worst <- vapply(
    rows,
    function(i) {
      smd <- i[x$statistic[i] == "smd" & !is.na(x$value[i])]
      if (length(smd) == 0) {
        return(NA_integer_)
      }
      smd[which.max(abs(x$value[smd]))]
    },
    integer(1)
  )
  energy <- vapply(
    rows,
    function(i) {
      en <- x$value[i][x$statistic[i] == "energy"]
      if (length(en) == 1) en else NA_real_
    },
    numeric(1)
  )
  first <- vapply(rows, `[`, integer(1), 1)
  out <- tibble::as_tibble(x[first, id_cols, drop = FALSE])
  out$max_abs_smd <- abs(x$value[worst])
  out$variable <- x$variable[worst]
  out$level <- x$level[worst]
  out$energy <- unname(energy)
  if (any(x$statistic == "variance_ratio")) {
    out$variance_ratio <- unname(vapply(
      rows,
      function(i) {
        vr <- x$value[i][x$statistic[i] == "variance_ratio"]
        vr <- vr[!is.na(vr)]
        if (length(vr) == 0) {
          return(NA_real_)
        }
        vr[which.max(abs(log(vr)))]
      },
      numeric(1)
    ))
  }
  out
}
