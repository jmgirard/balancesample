test_that("folds are rsample objects with the right sizes", {
  d <- sim_data()
  for (m in available_methods()) {
    set.seed(1)
    folds <- balanced_vfold_cv(d, v = 7, balance = c(x1, x2, cat), method = m)
    expect_s3_class(folds, c("balanced_vfold_cv", "rset"))
    expect_equal(nrow(folds), 7)
    expect_equal(folds$id, sprintf("Fold%d", 1:7))
    n_assess <- vapply(
      folds$splits,
      function(s) length(rsample::complement(s)),
      integer(1)
    )
    expect_true(all(n_assess %in% c(28L, 29L)), info = m)
    expect_true(assessed_once(folds, nrow(d)), info = m)
    expect_s3_class(rsample::analysis(folds$splits[[1]]), "data.frame")
  }
})

test_that("repeats give two ID columns and fresh partitions", {
  d <- sim_data()
  set.seed(1)
  folds <- balanced_vfold_cv(
    d,
    v = 5,
    repeats = 3,
    balance = c(x1, x2),
    method = "rerandomize"
  )
  expect_equal(nrow(folds), 15)
  expect_named(folds, c("splits", "id", "id2"))
  expect_match(pretty(folds), "repeated 3 times")
  expect_true(assessed_once(folds, nrow(d)))
  first <- rsample::complement(folds$splits[[1]])
  others <- lapply(folds$splits[6:10], rsample::complement)
  expect_false(any(vapply(others, setequal, logical(1), first)))
})

test_that("exact strata are balanced to within one row", {
  d <- sim_data()
  strata <- interaction(d$cat, d$sex, drop = TRUE)
  for (m in available_methods()) {
    set.seed(2)
    folds <- balanced_vfold_cv(
      d,
      v = 6,
      balance = c(x1, x2, x3),
      exact = c(cat, sex),
      method = m
    )
    expect_true(strata_exact(folds, strata), info = m)
    expect_true(assessed_once(folds, nrow(d)), info = m)
  }
})

test_that("exact strata follow the actual fold sizes", {
  d <- data.frame(s = factor(rep(letters[1:3], c(2, 5, 3))))
  for (seed in 1:20) {
    set.seed(seed)
    folds <- balanced_vfold_cv(d, v = 4, exact = s, method = "interaction")
    expect_true(shares_respected(folds, d$s))
    expect_true(strata_exact(folds, d$s))
  }
})

test_that("exact levels that contain the separator stay distinct", {
  d <- data.frame(
    a = rep(c("a:b", "a"), each = 2),
    b = rep(c("c", "b:c"), each = 2)
  )
  for (seed in 1:10) {
    set.seed(seed)
    folds <- balanced_vfold_cv(
      d,
      v = 2,
      exact = c(a, b),
      method = "interaction"
    )
    expect_equal(sort(rsample::assessment(folds$splits[[1]])$a), c("a", "a:b"))
  }
})

test_that("column names that match paste() arguments work", {
  for (nm in c("sep", "collapse")) {
    d <- data.frame(rep(c("a", "b"), each = 4), rep(1:2, 4))
    names(d) <- c(nm, "other")
    set.seed(1)
    folds <- balanced_vfold_cv(
      d,
      v = 2,
      exact = c(all_of(nm), other),
      method = "interaction"
    )
    expect_true(assessed_once(folds, nrow(d)), info = nm)
  }
})

test_that("exact alone works without balance columns", {
  d <- sim_data()
  set.seed(3)
  folds <- balanced_vfold_cv(d, v = 4, exact = cat)
  expect_true(strata_exact(folds, d$cat))
})

test_that("groups stay together", {
  d <- sim_data()
  for (m in available_methods()) {
    set.seed(4)
    folds <- balanced_vfold_cv(
      d,
      v = 5,
      balance = c(x1, x2),
      group = id,
      method = m
    )
    expect_true(groups_intact(folds, d$id), info = m)
    expect_true(assessed_once(folds, nrow(d)), info = m)
  }
})

test_that("groups and exact strata work together", {
  d <- sim_data()
  site <- factor(d$site)
  for (m in available_methods()) {
    set.seed(5)
    folds <- balanced_vfold_cv(
      d,
      v = 4,
      balance = x1,
      exact = site,
      group = id,
      method = m
    )
    expect_true(groups_intact(folds, d$id), info = m)
    # Exactness holds for groups (not rows) per stratum.
    n_assess <- vapply(
      folds$splits,
      function(s) {
        length(unique(d$id[rsample::complement(s)]))
      },
      integer(1)
    )
    n_groups <- vapply(
      folds$splits,
      function(s) {
        out <- rsample::complement(s)
        tapply(d$id[out], site[out], function(z) length(unique(z)))
      },
      numeric(2)
    )
    expect_true(
      all(apply(n_groups, 1, function(r) max(r) - min(r) <= 1)),
      info = m
    )
    m_groups <- tapply(d$id, site, function(z) length(unique(z)))
    share <- outer(as.vector(m_groups), n_assess) / sum(m_groups)
    expect_true(
      all(n_groups >= floor(share + 1e-9) & n_groups <= ceiling(share - 1e-9)),
      info = m
    )
  }
})

test_that("exact columns must be constant within groups", {
  d <- sim_data()
  expect_snapshot(
    balanced_vfold_cv(
      d,
      v = 4,
      balance = x1,
      exact = sex,
      group = id,
      method = "rerandomize"
    ),
    error = TRUE
  )
})

test_that("balanced folds beat random folds on the balance variables", {
  d <- sim_data(n = 300)
  set.seed(6)
  random <- mean(replicate(
    10,
    max_abs_smd(rsample::vfold_cv(d, v = 5), x1, x2, x3)
  ))
  for (m in setdiff(available_methods(), "interaction")) {
    set.seed(6)
    bal <- mean(replicate(
      5,
      max_abs_smd(
        balanced_vfold_cv(d, v = 5, balance = c(x1, x2, x3), method = m)
      )
    ))
    expect_lt(bal, random)
  }
})

test_that("rerandomize meets its threshold or warns", {
  d <- sim_data()
  set.seed(7)
  expect_no_warning(
    balanced_vfold_cv(d, v = 5, balance = c(x1, x2), method = "rerandomize")
  )
  expect_warning(
    balanced_vfold_cv(
      d,
      v = 5,
      balance = c(x1, x2),
      method = "rerandomize",
      options = list(p_accept = 1e-12, max_draws = 3)
    ),
    "threshold"
  )
})

test_that("method options are checked", {
  d <- sim_data()
  expect_snapshot(error = TRUE, {
    balanced_vfold_cv(
      d,
      balance = x1,
      method = "rerandomize",
      options = list(foo = 1)
    )
    balanced_vfold_cv(d, balance = x1, method = "nope")
    balanced_vfold_cv(
      d,
      balance = x1,
      method = "interaction",
      options = list(1)
    )
  })
})

test_that("anticlust options work and are checked", {
  d <- sim_data(n = 60)
  opts <- list(
    list(moments = 1),
    list(moments = 4, covariances = TRUE),
    list(k_neighbours = 5),
    list(local_maximum = FALSE),
    list(repetitions = 3)
  )
  for (o in opts) {
    set.seed(8)
    folds <- balanced_vfold_cv(
      d,
      v = 3,
      balance = c(x1, x2),
      exact = sex,
      options = o
    )
    expect_true(assessed_once(folds, nrow(d)), info = names(o))
    expect_true(strata_exact(folds, factor(d$sex)), info = names(o))
  }
  expect_snapshot(error = TRUE, {
    balanced_vfold_cv(d, balance = x1, options = list(moments = 5))
    balanced_vfold_cv(d, balance = x1, options = list(k_neighbours = 0.5))
    balanced_vfold_cv(d, balance = x1, options = list(objective = "diversity"))
  })
})

test_that("anticlust accepts square categorical balance data", {
  d <- data.frame(x = factor(letters[1:4]))
  set.seed(1)
  folds <- balanced_vfold_cv(d, v = 2, balance = x, options = list(moments = 1))
  expect_true(assessed_once(folds, 4))
})

test_that("kplus_features() builds each moment from the centered input", {
  x <- matrix(c(1, 2, 4, 8, 16, 32))
  cx <- x - mean(x)
  expected <- standardize_matrix(cbind(x, cx^2, cx^3, cx^4))
  expect_equal(kplus_features(x, 4, FALSE), expected)
  expect_equal(ncol(kplus_features(x, 3, FALSE)), 3)
})

test_that("rerandomize accepts close to p_accept with tied imbalance", {
  # Every split of four rows into two pairs has one of two imbalance values,
  # and the smaller one has probability 1/3.
  # With one draw allowed, the warning shows how often that draw was not
  # accepted: about 0.9 for an acceptance rate of 0.1, and 1/3 if every tie
  # with the threshold were accepted.
  x <- scale(matrix(c(0, 0, 1, 1)))
  set.seed(1)
  warned <- replicate(400, {
    w <- FALSE
    withCallingHandlers(
      method_rerandomize(x, c(2L, 2L), p_accept = 0.1, max_draws = 1),
      warning = function(cnd) {
        w <<- TRUE
        invokeRestart("muffleWarning")
      }
    )
    w
  })
  expect_gt(mean(warned), 0.8)
  expect_lt(mean(warned), 0.97)
})

test_that("local_maximum balances at least as well as one pass", {
  d <- sim_data(n = 120)
  x <- build_balance_matrix(d[c("x1", "x2", "x3")])
  f <- kplus_features(x, 2, FALSE)
  sizes <- sizes_from_props(nrow(d), rep(1 / 4, 4))
  for (seed in 1:5) {
    set.seed(seed)
    one <- method_anticlust(x, sizes, local_maximum = FALSE)
    set.seed(seed)
    many <- method_anticlust(x, sizes)
    expect_gte(
      anticlust::variance_objective(f, many),
      anticlust::variance_objective(f, one)
    )
  }
})

test_that("anticlust repetitions keep the exact counts", {
  d <- withr::with_seed(
    1,
    data.frame(x = stats::rnorm(19), s = factor(rep(1:3, c(5, 6, 8))))
  )
  set.seed(1)
  folds <- balanced_vfold_cv(
    d,
    4,
    balance = x,
    exact = s,
    options = list(repetitions = 20)
  )
  expect_true(shares_respected(folds, d$s))
  expect_true(strata_exact(folds, d$s))
})

test_that("anticlust handles factor-only and one-row-per-fold data", {
  set.seed(1)
  folds <- balanced_vfold_cv(data.frame(g = factor(1:6)), 3, balance = g)
  expect_true(assessed_once(folds, 6))
  folds <- balanced_vfold_cv(data.frame(x = 1:2), 2, balance = x)
  expect_true(assessed_once(folds, 2))
})

test_that("twinning falls back to random folds for small data", {
  skip_if_not_installed("twinning")
  set.seed(1)
  folds <- balanced_vfold_cv(
    data.frame(x = 1:5),
    3,
    balance = x,
    method = "twinning"
  )
  expect_true(assessed_once(folds, 5))
})

test_that("twinning checks strategy for two folds", {
  skip_if_not_installed("twinning")
  expect_error(
    balanced_vfold_cv(
      mtcars,
      2,
      balance = mpg,
      method = "twinning",
      options = list(strategy = 99)
    ),
    "strategy"
  )
  set.seed(1)
  folds <- balanced_vfold_cv(
    mtcars,
    2,
    balance = mpg,
    method = "twinning",
    options = list(strategy = 3)
  )
  expect_true(assessed_once(folds, nrow(mtcars)))
})

test_that("rerandomize calibrates its threshold under exact strata", {
  d <- data.frame(s = rep(0:1, each = 5))
  set.seed(1)
  expect_no_warning(
    balanced_initial_split(
      d,
      0.5,
      balance = s,
      exact = s,
      method = "rerandomize",
      options = list(max_draws = 5)
    )
  )
})

test_that("reshuffle_rset() works when balancesample is attached", {
  set.seed(1)
  folds <- balanced_vfold_cv(mtcars, 3, balance = mpg, method = "interaction")
  re <- rsample::reshuffle_rset(folds)
  expect_s3_class(re, "balanced_vfold_cv")
  expect_true(assessed_once(re, nrow(mtcars)))
})

test_that("twinning strategy 2 needs a power of two", {
  skip_if_not_installed("twinning")
  d <- sim_data()
  expect_snapshot(
    balanced_vfold_cv(
      d,
      v = 5,
      balance = x1,
      method = "twinning",
      options = list(strategy = 2)
    ),
    error = TRUE
  )
  set.seed(9)
  folds <- balanced_vfold_cv(
    d,
    v = 4,
    balance = x1,
    method = "twinning",
    options = list(strategy = 2)
  )
  expect_true(assessed_once(folds, nrow(d)))
})

test_that("bad inputs give clear errors", {
  d <- sim_data()
  expect_snapshot(error = TRUE, {
    balanced_vfold_cv(d)
    balanced_vfold_cv(d, v = 1, balance = x1)
    balanced_vfold_cv(d[1:3, ], v = 4, balance = x1)
    balanced_vfold_cv(d, balance = x1, exact = x2)
    balanced_vfold_cv(d, balance = c(x1, id), group = id)
    balanced_vfold_cv(d, balance = x1, group = c(id, sex))
    balanced_vfold_cv(d, v = 50, balance = x1, group = id)
  })
  d$x1[3] <- NA
  expect_snapshot(balanced_vfold_cv(d, balance = x1), error = TRUE)
  m <- data.frame(x = I(cbind(1:6, 7:12)))
  expect_snapshot(balanced_vfold_cv(m, 3, balance = x), error = TRUE)
})

test_that("strata and other extra arguments are rejected", {
  d <- sim_data()
  expect_snapshot(error = TRUE, {
    balanced_vfold_cv(d, balance = x1, strata = cat)
    balanced_initial_split(d, balance = x1, strata = cat)
    balanced_vfold_cv(d, balance = x1, folds = 5)
  })
})

test_that("results are reproducible with a seed", {
  d <- sim_data()
  for (m in available_methods()) {
    set.seed(10)
    a <- balanced_vfold_cv(d, v = 4, balance = c(x1, x2), method = m)
    set.seed(10)
    b <- balanced_vfold_cv(d, v = 4, balance = c(x1, x2), method = m)
    expect_identical(
      lapply(a$splits, rsample::complement),
      lapply(b$splits, rsample::complement),
      info = m
    )
  }
})

test_that("print shows the method and design", {
  d <- sim_data()
  set.seed(11)
  folds <- balanced_vfold_cv(
    d,
    v = 3,
    balance = x1,
    exact = site,
    group = id,
    method = "rerandomize"
  )
  expect_match(
    pretty(folds),
    "3-fold balanced cross-validation \\(rerandomize\\)"
  )
  expect_match(pretty(folds), "stratified by site, grouped by id")
  expect_output(print(folds), "balanced cross-validation")

  # A long header wraps between its parts.
  withr::local_options(width = 60)
  header <- utils::capture.output(print(folds))[1:2]
  expect_equal(
    header,
    c(
      "# 3-fold balanced cross-validation (rerandomize),",
      "# stratified by site, grouped by id"
    )
  )
  # A header that fits exactly stays on one line.
  full <- paste0("# ", pretty(folds))
  withr::local_options(width = nchar(full))
  expect_equal(utils::capture.output(print(folds))[1], full)
  withr::local_options(width = nchar(full) - 1)
  expect_equal(
    utils::capture.output(print(folds))[1:2],
    c(
      "# 3-fold balanced cross-validation (rerandomize), stratified by site,",
      "# grouped by id"
    )
  )
})

test_that("print wraps long lists of exact columns by display width", {
  d <- sim_data()
  d$sex2 <- d$sex
  d$site2 <- d$site
  # 16 wide characters, each two columns wide.
  names(d)[names(d) == "site2"] <- strrep("\u754c", 16)
  set.seed(12)
  folds <- balanced_vfold_cv(
    d,
    v = 3,
    balance = x1,
    exact = c(site, sex, sex2, tidyselect::last_col()),
    method = "rerandomize"
  )
  expect_match(pretty(folds), "stratified by site, sex, sex2, ")
  withr::local_options(width = 60)
  header <- grep("^# ", utils::capture.output(print(folds)), value = TRUE)
  header <- header[seq_len(length(header) - 1)]
  expect_true(all(nchar(header, type = "width") <= 60))
  expect_equal(header[2], "# stratified by site, sex, sex2,")
  expect_equal(
    paste(sub("^# ", "", header), collapse = " "),
    pretty(folds)
  )
})

test_that("bad data give clear errors", {
  expect_snapshot(error = TRUE, {
    balanced_vfold_cv(1:10, balance = x)
    balanced_vfold_cv(data.frame(x = 1), balance = x)
    balanced_vfold_cv(data.frame(x = c(1, Inf, 3, 4)), v = 2, balance = x)
    balanced_vfold_cv(
      data.frame(x = as.Date("2026-01-01") + 0:3),
      v = 2,
      balance = x
    )
  })
})

test_that("more anticlust options are checked", {
  d <- sim_data(n = 60)
  expect_snapshot(error = TRUE, {
    balanced_vfold_cv(d, balance = x1, options = list(covariances = "yes"))
    balanced_vfold_cv(d, balance = x1, options = list(local_maximum = NA))
    balanced_vfold_cv(d, balance = x1, options = list(repetitions = 0))
  })
})

test_that("rerandomize and interaction options are checked", {
  d <- sim_data(n = 60)
  expect_snapshot(error = TRUE, {
    balanced_vfold_cv(
      d,
      balance = x1,
      method = "rerandomize",
      options = list(p_accept = 0)
    )
    balanced_vfold_cv(
      d,
      balance = x1,
      method = "rerandomize",
      options = list(max_draws = 0)
    )
    balanced_vfold_cv(
      d,
      balance = x1,
      method = "interaction",
      options = list(n_bins = 1)
    )
    balanced_vfold_cv(
      d,
      balance = x1,
      method = "interaction",
      options = list(pool = 1)
    )
  })
})

test_that("logical balance columns and exact-only groups work", {
  d <- sim_data()
  d$flag <- d$x1 > 0
  set.seed(12)
  folds <- balanced_vfold_cv(d, v = 4, balance = c(flag, x2))
  expect_true(assessed_once(folds, nrow(d)))
  bc <- balance_check(folds)
  expect_equal(unique(bc$level[bc$variable %in% "flag"]), c("FALSE", "TRUE"))

  set.seed(13)
  folds <- balanced_vfold_cv(d, v = 4, exact = site, group = id)
  expect_true(groups_intact(folds, d$id))
  expect_true(assessed_once(folds, nrow(d)))
})

test_that("bad option values give the package's own errors", {
  classed <- structure(3, class = "classed_number")
  # A class can claim to be numeric; the checks must not ask it.
  local_mocked_s3_method("is.numeric", "classed_number", function(x) "yes")
  # Each case is a bad option and the start of its error message.
  must <- list(
    moments = "must be 1, 2, 3, or 4",
    covariances = "must be `TRUE` or `FALSE`",
    k_neighbours = "must be a whole number of at least 1, or",
    local_maximum = "must be `TRUE` or `FALSE`",
    repetitions = "must be a whole number of at least 1",
    strategy = "must be 1, 2, or 3",
    leaf_size = "must be a whole number of at least 1",
    p_accept = "must be a number above 0 and at most 1",
    max_draws = "must be a whole number of at least 1",
    n_bins = "must be a whole number of at least 2",
    pool = "must be a number from 0 up to"
  )
  bad <- list(
    anticlust = list(
      list(moments = NA_real_),
      list(moments = 2.5),
      list(moments = 5),
      list(moments = 0),
      list(moments = asS4(2)),
      list(moments = numeric(0)),
      list(covariances = matrix(TRUE)),
      list(covariances = 1),
      list(covariances = structure(TRUE, class = "classed_number")),
      list(k_neighbours = 0),
      list(k_neighbours = -Inf),
      list(k_neighbours = 1.5),
      list(k_neighbours = matrix(Inf)),
      list(k_neighbours = 1e20),
      list(k_neighbours = classed),
      list(local_maximum = NA),
      list(local_maximum = matrix(TRUE)),
      list(local_maximum = asS4(TRUE)),
      list(local_maximum = c(TRUE, FALSE)),
      list(repetitions = NA_real_),
      list(repetitions = 0),
      list(repetitions = 1.5)
    ),
    rerandomize = list(
      list(p_accept = NA_real_),
      list(p_accept = 2),
      list(p_accept = 0),
      list(p_accept = -0.1),
      list(p_accept = TRUE),
      list(p_accept = matrix(0.5)),
      list(p_accept = structure(0.5, class = "classed_number")),
      list(max_draws = NA_real_),
      list(max_draws = 2.5),
      list(max_draws = 2^52),
      list(max_draws = 0)
    ),
    interaction = list(
      list(n_bins = NA_real_),
      list(n_bins = 2.5),
      list(n_bins = matrix(2)),
      list(n_bins = 1),
      list(pool = NA_real_),
      list(pool = 1),
      list(pool = -0.1),
      list(pool = matrix(0.5))
    ),
    twinning = list(
      list(strategy = 1:2),
      list(strategy = 4),
      list(strategy = 0),
      list(leaf_size = NA_real_),
      list(leaf_size = 0),
      list(leaf_size = -1),
      list(leaf_size = 2.5),
      list(leaf_size = classed)
    )
  )
  # The twinning checks run without twinning, so test them directly too.
  for (o in bad$twinning) {
    nm <- names(o)
    err <- expect_error(
      rlang::exec(check_twinning_options, !!!o),
      class = "rlang_error"
    )
    expect_match(
      conditionMessage(err),
      paste0("`", nm, "` ", must[[nm]]),
      fixed = TRUE
    )
  }
  if (!rlang::is_installed("twinning")) {
    bad$twinning <- NULL
  }
  d <- data.frame(
    x = 1:12,
    const = 1,
    ex = rep(c("a", "b"), 6),
    g = rep(1:6, 2)
  )
  # With exact or constant columns only, the method doesn't run, but its
  # options are still checked.
  calls <- list(
    balance = function(m, o) {
      balanced_vfold_cv(d, 2, balance = x, method = m, options = o)
    },
    exact = function(m, o) {
      balanced_vfold_cv(d, 2, exact = ex, method = m, options = o)
    },
    constant = function(m, o) {
      balanced_vfold_cv(d, 2, balance = const, method = m, options = o)
    },
    group = function(m, o) {
      balanced_vfold_cv(d, 2, balance = x, group = g, method = m, options = o)
    },
    split = function(m, o) {
      balanced_initial_split(d, 1 / 2, balance = x, method = m, options = o)
    },
    validation = function(m, o) {
      balanced_initial_validation_split(
        d,
        c(0.5, 0.25),
        balance = x,
        method = m,
        options = o
      )
    }
  )
  for (m in names(bad)) {
    for (o in bad[[m]]) {
      nm <- names(o)
      for (cl in names(calls)) {
        # Two-way splits refuse `strategy` before the method is chosen.
        if (cl == "split" && nm == "strategy") {
          next
        }
        err <- expect_error(calls[[cl]](m, o), class = "rlang_error")
        expect_s3_class(err$parent, "rlang_error")
        expect_match(
          conditionMessage(err$parent),
          paste0("`", nm, "` ", must[[nm]]),
          fixed = TRUE,
          info = paste(m, nm, cl)
        )
      }
    }
  }
})

test_that("option names must not be missing", {
  expect_error(
    balanced_vfold_cv(
      mtcars,
      2,
      balance = mpg,
      options = stats::setNames(list(2), NA_character_)
    ),
    "unique names",
    class = "rlang_error"
  )
})
