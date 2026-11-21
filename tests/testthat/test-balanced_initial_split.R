test_that("initial splits have the right sizes and classes", {
  d <- sim_data()
  for (m in available_methods()) {
    set.seed(1)
    s <- balanced_initial_split(
      d,
      prop = 3 / 4,
      balance = c(x1, x2, cat),
      method = m
    )
    expect_s3_class(s, c("balanced_initial_split", "initial_split", "rsplit"))
    expect_equal(nrow(rsample::training(s)), 150, info = m)
    expect_equal(nrow(rsample::testing(s)), 50, info = m)
    expect_equal(attr(s, "prop"), 3 / 4)
    expect_equal(attr(s, "method"), m)
  }
})

test_that("initial splits stratify exactly", {
  d <- sim_data()
  strata <- interaction(d$cat, d$sex, drop = TRUE)
  share <- as.vector(table(strata)) * 0.8
  for (m in available_methods()) {
    set.seed(2)
    s <- balanced_initial_split(
      d,
      prop = 0.8,
      balance = c(x1, x2),
      exact = c(cat, sex),
      method = m
    )
    n_train <- tabulate(as.integer(strata)[s$in_id], nlevels(strata))
    expect_true(
      all(n_train >= floor(share) & n_train <= ceiling(share)),
      info = m
    )
  }
})

test_that("initial splits keep groups together", {
  d <- sim_data()
  for (m in available_methods()) {
    set.seed(3)
    s <- balanced_initial_split(
      d,
      prop = 0.75,
      balance = c(x1, x2),
      group = id,
      method = m
    )
    expect_length(intersect(d$id[s$in_id], d$id[rsample::complement(s)]), 0)
    # 38 groups * 0.75 = 28.5 groups for training
    expect_true(length(unique(d$id[s$in_id])) %in% 28:29, info = m)
  }
})

test_that("twinning needs prop = 1 - 1/r", {
  skip_if_not_installed("twinning")
  d <- sim_data()
  expect_snapshot(
    balanced_initial_split(d, prop = 0.7, balance = x1, method = "twinning"),
    error = TRUE
  )
  set.seed(4)
  s <- balanced_initial_split(d, prop = 0.25, balance = x1, method = "twinning")
  expect_equal(nrow(rsample::training(s)), 50)
})

test_that("two-way twinning splits reject `strategy`", {
  skip_if_not_installed("twinning")
  d <- sim_data()
  expect_snapshot(
    balanced_initial_split(
      d,
      balance = x1,
      method = "twinning",
      options = list(strategy = 2)
    ),
    error = TRUE
  )
})

test_that("twin_split_parts() depends on prop only", {
  expect_equal(twin_split_parts(3 / 4), c(3L, 1L))
  expect_equal(twin_split_parts(1 - 1 / 3), c(2L, 1L))
  expect_equal(twin_split_parts(0.9), c(9L, 1L))
  expect_equal(twin_split_parts(0.5), c(1L, 1L))
  expect_equal(twin_split_parts(0.25), c(1L, 3L))
  expect_equal(twin_split_parts(1 - 1 / 100000), c(99999, 1))
  expect_null(twin_split_parts(0.7))
  expect_null(twin_split_parts(0.67))
  expect_null(twin_split_parts(1e-310))
  expect_null(twin_split_parts(0.5001))
  # Typed decimals and proportions computed by subtraction pass.
  for (p in c(0.95, 0.99, 0.875, 0.05, 1 - 1 / 7, 1 / 7, 5 / 6, 1 - 0.9875)) {
    expect_false(is.null(twin_split_parts(p)))
  }
  expect_equal(twin_split_parts(1 - (1 - 1 / 100000)), c(1, 99999))
  # Values within the relative tolerance of 1e-8 count as the fraction.
  expect_equal(twin_split_parts(0.500000003), c(1, 1))
})

test_that("near misses give the same twinning split as the fraction", {
  skip_if_not_installed("twinning")
  d <- data.frame(x = stats::rnorm(199))
  sizes <- vapply(
    c(0.5, 0.500000003),
    function(p) {
      set.seed(4)
      s <- balanced_initial_split(d, prop = p, balance = x, method = "twinning")
      nrow(rsample::training(s))
    },
    numeric(1)
  )
  expect_equal(sizes[[1]], sizes[[2]])
})

test_that("twinning splits with a very large r fall back to random", {
  skip_if_not_installed("twinning")
  d <- sim_data()
  set.seed(5)
  s <- balanced_initial_split(
    d,
    prop = 1 / 2^32,
    balance = x1,
    method = "twinning"
  )
  expect_equal(nrow(rsample::training(s)), 1)
})

test_that("twinning accepts or rejects the same prop for any size", {
  skip_if_not_installed("twinning")
  d <- data.frame(x = 1:200)
  # With 10 rows, 0.7 gives sets of 7 and 3, which once passed as 1/r.
  for (n in c(10, 100, 200)) {
    expect_error(
      balanced_initial_split(
        d[seq_len(n), , drop = FALSE],
        prop = 0.7,
        balance = x,
        method = "twinning"
      ),
      "1/r of the data"
    )
    set.seed(n)
    s <- balanced_initial_split(
      d[seq_len(n), , drop = FALSE],
      prop = 2 / 3,
      balance = x,
      method = "twinning"
    )
    expect_equal(nrow(rsample::testing(s)), round(n / 3))
  }
})

test_that("bad prop values give clear errors", {
  d <- sim_data()
  expect_snapshot(error = TRUE, {
    balanced_initial_split(d, prop = 1, balance = x1)
    balanced_initial_split(d, prop = "a", balance = x1)
  })
})

test_that("small and extreme splits don't depend on the seed", {
  d <- data.frame(x = 1:6, g = rep(1:2, each = 3))
  for (seed in 1:10) {
    set.seed(seed)
    s <- balanced_initial_split(
      d,
      0.75,
      balance = x,
      group = g,
      method = "interaction"
    )
    expect_equal(length(s$in_id), 3)
  }
  set.seed(1)
  s <- balanced_initial_split(
    mtcars,
    0.999,
    balance = mpg,
    method = "rerandomize"
  )
  expect_equal(nrow(rsample::testing(s)), 1)
})

test_that("grouped initial splits need two groups", {
  d <- data.frame(x = 1:6, g = 1)
  expect_snapshot(
    balanced_initial_split(d, balance = x, group = g),
    error = TRUE
  )
})

test_that("print shows the method and the columns", {
  d <- sim_data()
  withr::local_options(width = 100)
  set.seed(1)
  s <- balanced_initial_split(d, balance = c(x1, x2), exact = site, group = id)
  out <- utils::capture.output(printed <- print(s))
  expect_equal(
    out[1:2],
    c(
      "# Balanced training/testing split (anticlust), balancing x1, x2, stratified by site, grouped by id",
      "<Training/Testing/Total>"
    )
  )
  expect_identical(printed, s)

  # The empty fallback of a calibration split has no header.
  empty <- rsample::make_splits(d, d[0, ])
  class(empty) <- class(s)
  expect_equal(
    utils::capture.output(print(empty))[1],
    "<Training/Testing/Total>"
  )
})
