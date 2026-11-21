# Column selection and preprocessing shared by the user-facing functions.

prepare_balance_data <- function(
  data,
  balance,
  exact,
  group,
  na = "error",
  call = rlang::caller_env()
) {
  na <- rlang::arg_match(na, c("error", "include"), error_call = call)
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "{.arg data} must be a data frame, not {.obj_type_friendly {data}}.",
      call = call
    )
  }
  if (nrow(data) < 2) {
    cli::cli_abort("{.arg data} must have at least two rows.", call = call)
  }

  bal_pos <- tidyselect::eval_select(
    balance,
    data,
    allow_rename = FALSE,
    error_call = call
  )
  exa_pos <- tidyselect::eval_select(
    exact,
    data,
    allow_rename = FALSE,
    error_call = call
  )
  grp_pos <- tidyselect::eval_select(
    group,
    data,
    allow_rename = FALSE,
    error_call = call
  )

  if (length(bal_pos) == 0 && length(exa_pos) == 0) {
    cli::cli_abort(
      c(
        "Select at least one column with {.arg balance} or {.arg exact}.",
        i = "For example, {.code balance = c(age, income)}."
      ),
      call = call
    )
  }
  if (length(grp_pos) > 1) {
    cli::cli_abort(
      "{.arg group} must select one column, not {length(grp_pos)}.",
      call = call
    )
  }
  if (length(grp_pos) == 1 && any(c(bal_pos, exa_pos) == grp_pos)) {
    cli::cli_abort(
      "The {.arg group} column {.field {names(grp_pos)}} can't also be in {.arg balance} or {.arg exact}.",
      call = call
    )
  }

  used <- unique(c(bal_pos, exa_pos, grp_pos))
  not_vector <- names(data)[used][vapply(
    data[used],
    function(v) !is.null(dim(v)) || is.list(v),
    logical(1)
  )]
  if (length(not_vector) > 0) {
    is_lt <- vapply(data[not_vector], inherits, logical(1), "POSIXlt")
    cli::cli_abort(
      c(
        "Selected columns must be plain vectors.",
        x = "{.field {not_vector}} {?is a/are} matrix, data frame, or list column{?s}.",
        i = if (any(is_lt)) {
          "Convert {.cls POSIXlt} date-times with {.fn as.POSIXct}."
        }
      ),
      call = call
    )
  }
  check_missing(data, bal_pos, exa_pos, grp_pos, na, call = call)

  for (nm in names(exa_pos)) {
    col <- data[[nm]]
    if (is.numeric(col) && any(is.infinite(col))) {
      cli::cli_abort("{.field {nm}} has infinite values.", call = call)
    }
    whole <- is.numeric(col) && all(col == round(col), na.rm = TRUE)
    if (!(is.factor(col) || is.character(col) || is.logical(col) || whole)) {
      cli::cli_abort(
        c(
          "{.arg exact} columns must be categorical.",
          x = "{.field {nm}} is {.obj_type_friendly {col}}.",
          i = "Convert it to a factor, or move it to {.arg balance}."
        ),
        call = call
      )
    }
  }

  exact_fct <- NULL
  if (length(exa_pos) > 0) {
    exact_fct <- combine_levels(data[exa_pos])
  }
  group_fct <- NULL
  if (length(grp_pos) == 1) {
    group_fct <- factor(data[[grp_pos]])
  }
  if (!is.null(exact_fct) && !is.null(group_fct)) {
    check_exact_in_groups(exact_fct, group_fct, data[exa_pos], call = call)
  }

  list(
    x = build_balance_matrix(data[bal_pos], call = call),
    exact = exact_fct,
    group = group_fct,
    balance_names = names(bal_pos),
    exact_names = names(exa_pos),
    exact_numeric = vapply(data[exa_pos], is.numeric, logical(1)),
    group_name = names(grp_pos),
    na = na
  )
}

# With `na = "error"`, no selected column may have missing values. With
# `na = "include"`, `balance` and `exact` columns may, as long as each has an
# observed value. `group` never may: rows with a missing ID may or may not
# belong together, and only the user knows which.
check_missing <- function(data, bal_pos, exa_pos, grp_pos, na, call) {
  used <- unique(c(bal_pos, exa_pos, grp_pos))
  n_na <- vapply(data[used], function(v) sum(is.na(v)), integer(1))
  if (!any(n_na > 0)) {
    return(invisible())
  }
  group_hint <- if (length(grp_pos) == 1 && n_na[[names(grp_pos)]] > 0) {
    "{.arg group} can't have missing values. Give the rows with a missing ID their own IDs, or one shared ID, depending on whether they can belong together."
  }
  if (na == "error") {
    n <- n_na[n_na > 0]
    has_na <- paste0(names(n), " (", n, ifelse(n == 1, " row)", " rows)"))
    cli::cli_abort(
      c(
        "Selected columns must not have missing values.",
        x = "Missing values in {has_na}.",
        i = if (any(n[setdiff(names(n), names(grp_pos))] > 0)) {
          "Use {.code na = \"include\"} to balance on which values are missing, or drop those rows. See {.code vignette(\"balancesample\")}."
        },
        i = group_hint
      ),
      call = call
    )
  }
  if (!is.null(group_hint)) {
    n <- n_na[[names(grp_pos)]]
    cli::cli_abort(
      c(
        "{.arg group} must not have missing values.",
        x = "{.field {names(grp_pos)}} is missing in {n} row{?s}.",
        i = "Give the rows with a missing ID their own IDs, or one shared ID, depending on whether they can belong together."
      ),
      call = call
    )
  }
  empty <- names(n_na)[n_na == nrow(data)]
  if (length(empty) > 0) {
    cli::cli_abort(
      c(
        "{.arg balance} and {.arg exact} columns must have at least one observed value.",
        x = "{.field {empty}} {?is/are} missing in every row."
      ),
      call = call
    )
  }
  invisible()
}

# `exact` columns must be constant within each group, so that each group has
# one stratum. A missing value counts as its own level.
check_exact_in_groups <- function(exact, group, exact_cols, call) {
  g <- as.integer(group)
  n_lvl <- vapply(
    split(as.integer(exact), g),
    function(z) length(unique(z)),
    integer(1)
  )
  if (all(n_lvl <= 1)) {
    return(invisible())
  }
  # Does a group that fails have both missing and observed values in some
  # column?
  mixed_na <- any(vapply(
    exact_cols,
    function(v) {
      some <- tapply(is.na(v), g, function(z) any(z) && !all(z))
      any(some[n_lvl > 1])
    },
    logical(1)
  ))
  cli::cli_abort(
    c(
      "{.arg exact} columns must be constant within each group.",
      x = "{sum(n_lvl > 1)} group{?s} contain{?s/} more than one {.arg exact} level.",
      i = if (mixed_na) {
        "A missing value counts as its own level. Fill it in where you know the group's value."
      }
    ),
    call = call
  )
}

# One factor level per observed combination of the columns in `cols`. The key
# is built from integer codes, so level labels that contain the separator
# can't make two combinations collide. A missing value gets a code of its own.
combine_levels <- function(cols) {
  codes <- lapply(cols, function(v) {
    as.integer(addNA(as.factor(v), ifany = TRUE))
  })
  factor(do.call(paste, c(unname(codes), sep = "_")))
}

# Numeric matrix with one standardized column per numeric variable and one
# standardized indicator column per level of each categorical variable.
# Missing categorical values get an indicator column of their own. Missing
# numeric values get the mean of the observed values, and an indicator
# column marks them. Constant columns are dropped.
build_balance_matrix <- function(df, call = rlang::caller_env()) {
  n <- nrow(df)
  if (ncol(df) == 0) {
    return(matrix(numeric(0), nrow = n, ncol = 0))
  }
  cols <- lapply(names(df), function(nm) {
    v <- df[[nm]]
    if (is.logical(v)) {
      v <- as.integer(v)
    }
    if (is.character(v)) {
      v <- factor(v)
    }
    if (is.factor(v)) {
      v <- droplevels(v)
      codes <- as.integer(v)
      lv <- levels(v)
      if (anyNA(codes)) {
        codes[is.na(codes)] <- length(lv) + 1L
        lv <- c(lv, "(missing)")
      }
      m <- outer(codes, seq_along(lv), "==") * 1
      colnames(m) <- paste0(nm, "_", lv)
      return(m)
    }
    v <- as_balance_numeric(v)
    if (is.numeric(v)) {
      if (any(is.infinite(v))) {
        cli::cli_abort("{.field {nm}} has infinite values.", call = call)
      }
      v <- as.double(v)
      miss <- is.na(v)
      if (all(miss)) {
        return(matrix(numeric(0), nrow = n, ncol = 0))
      }
      # Measure from an observed value, so that a large origin (such as a
      # date-time in seconds) doesn't swamp small differences.
      v <- v - v[!miss][1]
      if (!any(miss)) {
        return(matrix(v, ncol = 1, dimnames = list(NULL, nm)))
      }
      v[miss] <- mean(v[!miss])
      m <- cbind(v, miss * 1)
      colnames(m) <- c(nm, paste0(nm, "_(missing)"))
      return(m)
    }
    cli::cli_abort(
      c(
        "{.arg balance} columns must be numeric, logical, character, factor, Date, POSIXct, or difftime.",
        x = "{.field {nm}} is {.obj_type_friendly {v}}."
      ),
      call = call
    )
  })
  standardize_matrix(do.call(cbind, cols))
}

# Dates count in days, date-times in seconds, and durations in seconds, so
# that durations in different units give the same values. Other vectors are
# returned unchanged.
as_balance_numeric <- function(v) {
  if (inherits(v, "difftime")) {
    return(as.numeric(v, units = "secs"))
  }
  if (inherits(v, c("Date", "POSIXct"))) {
    return(as.numeric(v))
  }
  v
}

# Columns whose standard deviation is within rounding error of zero are
# dropped. By default, the tolerance is relative to the largest distance from
# the column's first value, so it depends on neither the units nor the
# origin. `size` gives the scale of the values that the columns were
# computed from instead.
standardize_matrix <- function(x, size = NULL) {
  if (ncol(x) == 0) {
    return(x)
  }
  if (is.null(size)) {
    size <- apply(x, 2, function(col) max(abs(col - col[1])))
  }
  sds <- apply(x, 2, stats::sd)
  keep <- is.finite(sds) & sds > near_zero(size)
  x <- x[, keep, drop = FALSE]
  if (ncol(x) == 0) {
    return(x)
  }
  centers <- colMeans(x)
  x <- sweep(x, 2, centers, "-")
  sweep(x, 2, sds[keep], "/")
}

# Collapse rows to one row per group, assign groups, and expand back to rows.
# Each group row holds the group means of the standardized columns plus the log
# group size, so large and small groups are also spread over the sets.
assign_grouped <- function(
  x,
  exact,
  group,
  props,
  assign_fun
) {
  g <- as.integer(group)
  n_grp <- nlevels(group)
  counts <- tabulate(g, n_grp)

  if (ncol(x) > 0) {
    xg <- rowsum(x, g, reorder = TRUE) / counts
  } else {
    xg <- matrix(numeric(0), nrow = n_grp, ncol = 0)
  }
  # Rounding error in a group mean is relative to the row values, not to the
  # means, which can all be equal up to rounding.
  size <- if (ncol(x) > 0) apply(abs(x), 2, max) else numeric(0)
  if (ncol(x) > 0 || !is.null(exact)) {
    xg <- cbind(xg, .group_size = log(counts))
    size <- c(size, max(log(counts)))
  }
  xg <- standardize_matrix(xg, size)

  exact_g <- NULL
  if (!is.null(exact)) {
    first <- match(seq_len(n_grp), g)
    exact_g <- droplevels(exact[first])
  }

  sizes_g <- sizes_from_props(n_grp, props)
  out_g <- assign_fun(xg, sizes_g, exact_g)
  out_g[g]
}

# Rounding error for values up to `size` in absolute value.
near_zero <- function(size) {
  64 * .Machine$double.eps * size
}

# Warn when more than half of the `exact` strata are missing from at least
# one set of the assignment `labels`. Strata are counted in groups with
# `group`. A stratum with at least n / (size of the smallest set) units has a
# share of one or more in every set, so the exact counts put it in every set.
# Smaller strata can be in every set too, so the message gives this count as
# enough, not as needed.
warn_sparse_strata <- function(prep, labels, k, call = rlang::caller_env()) {
  if (is.null(prep$exact)) {
    return(invisible())
  }
  unit <- "rows"
  strata <- as.integer(droplevels(prep$exact))
  if (!is.null(prep$group)) {
    unit <- "groups"
    first <- !duplicated(prep$group)
    strata <- strata[first]
    labels <- labels[first]
  }
  n_sets <- tapply(labels, strata, function(z) length(unique(z)))
  short <- sum(n_sets < k)
  if (short <= length(n_sets) / 2) {
    return(invisible())
  }
  needed <- ceiling(length(labels) / min(tabulate(labels, k)) - 1e-9)
  numeric_exact <- prep$exact_names[prep$exact_numeric]
  cli::cli_warn(
    c(
      "{short} of {length(n_sets)} {.arg exact} strata are missing from at least one set.",
      i = "A stratum with at least {needed} {unit} is in every set.",
      i = if (length(numeric_exact) > 0) {
        "Whole-number columns in {.arg exact} are categories, and each value is a level: {.field {numeric_exact}}. Put continuous columns in {.arg balance}, or use {.fn cut}."
      }
    ),
    class = "balancesample_sparse_strata",
    call = call
  )
  invisible()
}
