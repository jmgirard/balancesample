test_that("validation splits have the right sizes and classes", {
  d <- sim_data()
  for (m in available_methods()) {
    set.seed(1)
    s <- balanced_initial_validation_split(
      d,
      prop = c(0.6, 0.2),
      balance = c(x1, x2, cat),
      method = m
    )
    expect_s3_class(
      s,
      c(
        "balanced_initial_validation_split",
        "initial_validation_split",
        "three_way_split"
      )
    )
    expect_equal(nrow(rsample::training(s)), 120, info = m)
    expect_equal(nrow(rsample::validation(s)), 40, info = m)
    expect_equal(nrow(rsample::testing(s)), 40, info = m)
    expect_length(intersect(s$train_id, s$val_id), 0)
    expect_equal(attr(s, "val_att")$prop, c(0.6, 0.2))
    expect_equal(attr(s, "val_att")$method, m)
  }
})

test_that("validation splits stratify exactly", {
  d <- sim_data()
  strata <- interaction(d$cat, d$sex, drop = TRUE)
  m <- as.vector(table(strata))
  for (meth in available_methods()) {
    set.seed(2)
    s <- balanced_initial_validation_split(
      d,
      prop = c(0.5, 0.25),
      balance = c(x1, x2),
      exact = c(cat, sex),
      method = meth
    )
    test_id <- setdiff(seq_len(nrow(d)), c(s$train_id, s$val_id))
    for (ids in list(s$train_id, s$val_id, test_id)) {
      cnt <- tabulate(as.integer(strata)[ids], nlevels(strata))
      share <- m * length(ids) / nrow(d)
      expect_true(
        all(cnt >= floor(share + 1e-9) & cnt <= ceiling(share - 1e-9)),
        info = meth
      )
    }
  }
})

test_that("validation splits keep groups together", {
  d <- sim_data()
  for (m in available_methods()) {
    set.seed(3)
    s <- balanced_initial_validation_split(
      d,
      prop = c(0.5, 0.25),
      balance = c(x1, x2),
      group = id,
      method = m
    )
    test_id <- setdiff(seq_len(nrow(d)), c(s$train_id, s$val_id))
    sets <- list(d$id[s$train_id], d$id[s$val_id], d$id[test_id])
    expect_length(intersect(sets[[1]], sets[[2]]), 0)
    expect_length(intersect(sets[[1]], sets[[3]]), 0)
    expect_length(intersect(sets[[2]], sets[[3]]), 0)
    # 38 groups * 0.5 = 19 groups for training
    expect_equal(length(unique(sets[[1]])), 19, info = m)
  }
})

test_that("validation splits balance better than rsample", {
  d <- sim_data()
  set.seed(4)
  s <- balanced_initial_validation_split(d, balance = c(x1, x2, x3))
  set.seed(4)
  r <- rsample::initial_validation_split(d)
  bc_s <- balance_check(s)
  bc_r <- balance_check(r, x1, x2, x3)
  expect_equal(unique(bc_s$id), c("Validation", "Testing"))
  expect_lt(
    max(abs(bc_s$value[bc_s$statistic == "smd"])),
    max(abs(bc_r$value[bc_r$statistic == "smd"]))
  )
})

test_that("twinning validation splits use whole parts", {
  skip_if_not_installed("twinning")
  d <- sim_data()
  set.seed(5)
  s <- balanced_initial_validation_split(
    d,
    prop = c(0.5, 0.25),
    balance = c(x1, x2),
    method = "twinning",
    options = list(strategy = 2)
  )
  expect_equal(c(length(s$train_id), length(s$val_id)), c(100, 50))
  expect_snapshot(error = TRUE, {
    balanced_initial_validation_split(
      d,
      prop = c(0.62, 0.2),
      balance = x1,
      method = "twinning"
    )
    balanced_initial_validation_split(
      d,
      prop = c(0.6, 0.2),
      balance = x1,
      method = "twinning",
      options = list(strategy = 2)
    )
  })
})

test_that("twin_parts() depends on the proportions only", {
  expect_equal(twin_parts(c(0.6, 0.2, 0.2)), c(3L, 1L, 1L))
  expect_equal(twin_parts(c(0.7, 0.15, 0.15)), c(14L, 3L, 3L))
  expect_equal(twin_parts(c(0.125, 0.125, 0.75)), c(1L, 1L, 6L))
  expect_null(twin_parts(c(0.62, 0.2, 0.18)))
})

test_that("twinning validation splits accept or reject the same proportions for any size", {
  skip_if_not_installed("twinning")
  d <- data.frame(x = 1:80)
  set.seed(8)
  s <- balanced_initial_validation_split(
    d[1:16, , drop = FALSE],
    prop = c(0.125, 0.125),
    balance = x,
    method = "twinning",
    options = list(strategy = 2)
  )
  expect_equal(c(length(s$train_id), length(s$val_id)), c(2, 2))
  for (n in c(50, 80)) {
    expect_error(
      balanced_initial_validation_split(
        d[seq_len(n), , drop = FALSE],
        prop = c(0.62, 0.2),
        balance = x,
        method = "twinning"
      ),
      "multiple of 1/m"
    )
  }
})

test_that("bad validation split inputs give clear errors", {
  d <- sim_data()
  expect_snapshot(error = TRUE, {
    balanced_initial_validation_split(d, prop = 0.6, balance = x1)
    balanced_initial_validation_split(d, prop = c(0.6, 0.4), balance = x1)
    balanced_initial_validation_split(d, balance = x1, strata = cat)
    balanced_initial_validation_split(
      d[1:4, ],
      balance = x1,
      group = site
    )
  })
})

test_that("validation sets work with tune", {
  # tune 1.1.2 is the first version whose last_fit() takes three-way splits.
  skip_if_not_installed("tune", "1.1.2")
  skip_if_not_installed("parsnip")

  set.seed(6)
  s <- balanced_initial_validation_split(
    mtcars,
    balance = c(mpg, hp, wt),
    method = "rerandomize"
  )
  vs <- rsample::validation_set(s)
  expect_equal(attr(vs, "balance"), c("mpg", "hp", "wt"))
  expect_equal(balance_check(vs)$id[1], "validation")

  res <- tune::fit_resamples(
    parsnip::linear_reg(),
    mpg ~ hp + wt,
    resamples = vs
  )
  expect_true(all(tune::collect_metrics(res)$n == 1))

  res <- tune::last_fit(parsnip::linear_reg(), mpg ~ hp + wt, split = s)
  expect_equal(nrow(tune::collect_predictions(res)), nrow(rsample::testing(s)))
})

test_that("calibration splits accept classed numeric proportions", {
  set.seed(11)
  s <- balanced_initial_validation_split(
    mtcars,
    prop = I(c(0.5, 0.35)),
    balance = hp,
    method = "rerandomize"
  )
  expect_no_warning(
    cal <- rsample::internal_calibration_split(s, split_args = list())
  )
  expect_equal(nrow(rsample::calibration(cal)), 7)
})

test_that("calibration splits of a validation split balance the training set", {
  d <- sim_data()
  set.seed(7)
  s <- balanced_initial_validation_split(
    d,
    prop = c(0.5, 0.25),
    balance = c(x1, x2),
    group = id
  )
  cal <- rsample::internal_calibration_split(s, split_args = list())
  expect_s3_class(cal, "balanced_initial_validation_split_cal")
  expect_s3_class(cal, "internal_calibration_split")
  fit_set <- rsample::analysis(cal)
  cal_set <- rsample::calibration(cal)
  expect_equal(nrow(fit_set) + nrow(cal_set), length(s$train_id))
  expect_length(intersect(fit_set$id, cal_set$id), 0)
  # 19 training groups, 1/3 of them for calibration
  expect_equal(length(unique(cal_set$id)), 19 / 3, tolerance = 1)
})

test_that("twinning calibration splits use the nearest share twinning can make", {
  expect_equal(nearest_twin_share(14, 3), 5 / 6)
  # 1/r is not evenly spaced: 7/17 is nearer 1/3 than 1/2.
  expect_equal(nearest_twin_share(10, 7), 2 / 3)
  expect_equal(nearest_twin_share(3, 1), 3 / 4)
  expect_equal(nearest_twin_share(1, 1), 0.5)
  expect_equal(nearest_twin_share(1, 4), 0.2)
  # 5/12 is exactly between 1/2 and 1/3. The tie goes to the larger r.
  expect_equal(nearest_twin_share(7, 5), 2 / 3)
  expect_equal(nearest_twin_share(5, 7), 1 / 3)
  # Large values don't overflow.
  expect_equal(nearest_twin_share(1e308, 6e307), 2 / 3)
  # 7/24 is exactly between 1/3 and 1/4, and stays a tie when scaled.
  expect_equal(nearest_twin_share(17 * 2^500, 7 * 2^500), 3 / 4)
  expect_equal(
    three_way_fit_share(list(
      prop = c(1500000000L, 1000000000L),
      method = "twinning"
    )),
    2 / 3
  )

  skip_if_not_installed("twinning")
  d <- sim_data()
  set.seed(9)
  s <- balanced_initial_validation_split(
    d,
    prop = c(0.7, 0.15),
    balance = c(x1, x2),
    method = "twinning"
  )
  # 140 training rows; 14/17 can't be twinned, so 1/6 go to calibration.
  expect_no_warning(
    cal <- rsample::internal_calibration_split(s, split_args = list())
  )
  expect_equal(nrow(rsample::calibration(cal)), round(140 / 6))

  # A share exactly between two twinning shares rounds as plain arithmetic
  # does. 7/12 is between 1/2 and 2/3.
  set.seed(10)
  s7 <- balanced_initial_validation_split(
    sim_data(320),
    prop = c(7, 5) / 16,
    balance = x1,
    method = "twinning"
  )
  cal <- rsample::internal_calibration_split(s7, split_args = list())
  expect_equal(nrow(rsample::calibration(cal)), round(140 / 3))

  # A near miss gives the same calibration split as the fraction it matches.
  cal_rows <- vapply(
    list(c(7, 5) / 13, c(7, 5) / 13 + c(1, -1) * 1e-10),
    function(p) {
      set.seed(4)
      s13 <- balanced_initial_validation_split(
        sim_data(199),
        prop = p,
        balance = x1,
        method = "twinning"
      )
      set.seed(5)
      cal <- rsample::internal_calibration_split(s13, split_args = list())
      nrow(rsample::calibration(cal))
    },
    numeric(1)
  )
  expect_equal(cal_rows[[1]], cal_rows[[2]])

  # Very large values don't overflow.
  cal <- rsample::internal_calibration_split(
    s,
    split_args = list(prop = rep(.Machine$double.xmax, 2))
  )
  expect_equal(nrow(rsample::calibration(cal)), 70)

  # A named method string is still twinning.
  attr(s, "val_att")$method <- c(chosen = "twinning")
  cal <- rsample::internal_calibration_split(s, split_args = list())
  expect_equal(nrow(rsample::calibration(cal)), round(140 / 6))

  # Bad shares still reach the fallback.
  expect_warning(
    cal <- rsample::internal_calibration_split(
      s,
      split_args = list(prop = c(NA_real_, 1))
    ),
    "empty calibration set"
  )
  expect_equal(nrow(rsample::calibration(cal)), 0)
  bad <- list(
    list(prop = c("0.5", "0.35")),
    list(prop = numeric(0)),
    list(prop = c(-1, -1)),
    list(prop = c(0.3, 0.2, 0.1)),
    list(options = identity),
    list(method = new.env())
  )
  for (args in bad) {
    expect_warning(
      cal <- rsample::internal_calibration_split(s, split_args = args),
      "empty calibration set"
    )
    expect_equal(nrow(rsample::calibration(cal)), 0)
  }

  # Errors from classed prop and options values reach the fallback too.
  local_mocked_s3_method("Ops", "bad_prop", function(e1, e2) stop("no math"))
  local_mocked_s3_method("$<-", "bad_options", function(x, name, value) {
    stop("read only")
  })
  classed <- list(
    "no math" = list(prop = structure(c(0.5, 0.35), class = "bad_prop")),
    "read only" = list(options = structure(list(), class = "bad_options"))
  )
  for (cause in names(classed)) {
    w <- expect_warning(
      cal <- rsample::internal_calibration_split(s, classed[[cause]]),
      "empty calibration set"
    )
    expect_match(conditionMessage(w$parent), cause)
    expect_equal(nrow(rsample::calibration(cal)), 0)
  }
})
