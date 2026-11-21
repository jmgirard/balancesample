# Registry of balancing methods.
#
# Each method is a function `(x_scaled, sizes, exact, group, ...)` that
# returns an integer vector of set labels (1 to length(sizes)) with one entry
# per row of `x_scaled`. When `group` is NULL, the label counts must equal
# `sizes`. Methods with `native_group = FALSE` never see `group`: the driver
# collapses groups to one row each before it calls them. `check` takes the
# options and errors on bad values. The driver calls it when it doesn't run
# the method (when there are no balance columns), so bad options are caught
# either way.

balance_method_registry <- new.env(parent = emptyenv())

register_balance_method <- function(
  name,
  fun,
  options = character(),
  native_group = FALSE,
  package = NULL,
  check = function(...) invisible()
) {
  balance_method_registry[[name]] <- list(
    name = name,
    fun = fun,
    options = options,
    native_group = native_group,
    package = package,
    check = check
  )
  invisible(name)
}

balance_method_names <- function() {
  sort(ls(balance_method_registry))
}

get_balance_method <- function(method, call = rlang::caller_env()) {
  if (!rlang::is_string(method)) {
    cli::cli_abort("{.arg method} must be a single string.", call = call)
  }
  method <- rlang::arg_match(
    method,
    balance_method_names(),
    error_arg = "method",
    error_call = call
  )
  entry <- balance_method_registry[[method]]
  if (!is.null(entry$package)) {
    rlang::check_installed(
      entry$package,
      reason = sprintf("to use `method = \"%s\"`.", method),
      call = call
    )
  }
  entry
}

check_method_options <- function(entry, options, call = rlang::caller_env()) {
  if (!is.list(options)) {
    cli::cli_abort(
      "{.arg options} must be a list, not {.obj_type_friendly {options}}.",
      call = call
    )
  }
  if (length(options) == 0) {
    return(invisible(options))
  }
  nms <- names(options)
  if (is.null(nms) || anyNA(nms) || any(nms == "") || anyDuplicated(nms)) {
    cli::cli_abort(
      "All elements of {.arg options} must have unique names.",
      call = call
    )
  }
  bad <- setdiff(nms, entry$options)
  if (length(bad) > 0) {
    valid <- entry$options
    cli::cli_abort(
      c(
        "Unknown {cli::qty(length(bad))}option{?s} for method {.val {entry$name}}: {.arg {bad}}.",
        i = if (length(valid)) {
          "Valid options: {.arg {valid}}."
        } else {
          "This method has no options."
        }
      ),
      call = call
    )
  }
  invisible(options)
}

# Run one method once and return row-level labels.
assign_balanced <- function(
  prep,
  props,
  entry,
  options,
  call = rlang::caller_env()
) {
  run <- function(x, sizes, exact, group = NULL) {
    rlang::try_fetch(
      if (ncol(x) == 0) {
        rlang::exec(entry$check, !!!options)
        stratified_random(sizes, exact)
      } else {
        rlang::exec(entry$fun, x, sizes, exact, group, !!!options)
      },
      error = function(cnd) {
        cli::cli_abort(
          "Method {.val {entry$name}} failed.",
          parent = cnd,
          call = call
        )
      }
    )
  }

  k <- length(props)
  if (!is.null(prep$group) && !entry$native_group) {
    out <- assign_grouped(
      prep$x,
      prep$exact,
      prep$group,
      props,
      run
    )
  } else {
    sizes <- sizes_from_props(nrow(prep$x), props)
    out <- run(prep$x, sizes, prep$exact, prep$group)
    if (is.null(prep$group) && !identical(tabulate(out, k), sizes)) {
      cli::cli_abort(
        "Internal error: method {.val {entry$name}} returned the wrong set sizes.",
        .internal = TRUE
      )
    }
  }

  out <- as.integer(out)
  if (length(out) != nrow(prep$x) || anyNA(out) || any(out < 1L | out > k)) {
    cli::cli_abort(
      "Internal error: method {.val {entry$name}} returned invalid labels.",
      .internal = TRUE
    )
  }
  if (any(tabulate(out, k) == 0)) {
    cli::cli_abort(
      "At least one set is empty. Use fewer sets or more data.",
      call = call
    )
  }
  out
}
