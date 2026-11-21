.onLoad <- function(libname, pkgname) {
  register_balance_method(
    "anticlust",
    method_anticlust,
    options = c(
      "moments",
      "covariances",
      "k_neighbours",
      "local_maximum",
      "repetitions"
    ),
    check = check_anticlust_options
  )
  register_balance_method(
    "twinning",
    method_twinning,
    options = c("strategy", "leaf_size"),
    package = "twinning",
    check = check_twinning_options
  )
  register_balance_method(
    "rerandomize",
    method_rerandomize,
    options = c("p_accept", "max_draws"),
    check = check_rerandomize_options
  )
  register_balance_method(
    "interaction",
    method_interaction,
    options = c("n_bins", "pool"),
    check = check_interaction_options
  )
  invisible()
}
