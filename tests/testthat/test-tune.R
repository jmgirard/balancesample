test_that("balanced folds work with tune::fit_resamples()", {
  # The floor in DESCRIPTION; much older tune can't take a model spec here.
  skip_if_not_installed("tune", "2.0.0")
  skip_if_not_installed("parsnip")

  set.seed(1)
  folds <- balanced_vfold_cv(
    mtcars,
    v = 4,
    balance = c(mpg, hp, wt),
    exact = am,
    method = "rerandomize"
  )
  res <- tune::fit_resamples(
    parsnip::linear_reg(),
    mpg ~ hp + wt,
    resamples = folds
  )
  expect_s3_class(res, "resample_results")
  metrics <- tune::collect_metrics(res)
  expect_equal(sort(metrics$.metric), c("rmse", "rsq"))
  expect_true(all(metrics$n == 4))
})

test_that("balanced initial splits work with tune::last_fit()", {
  skip_if_not_installed("tune", "2.0.0")
  skip_if_not_installed("parsnip")

  set.seed(2)
  split <- balanced_initial_split(
    mtcars,
    prop = 3 / 4,
    balance = c(mpg, hp, wt),
    method = "rerandomize"
  )
  res <- tune::last_fit(parsnip::linear_reg(), mpg ~ hp + wt, split = split)
  expect_s3_class(res, "last_fit")
  expect_equal(nrow(tune::collect_predictions(res)), 8)
})

test_that("calibration splits of folds are balanced initial splits", {
  d <- sim_data()
  set.seed(3)
  folds <- balanced_vfold_cv(
    d,
    v = 5,
    balance = c(x1, x2),
    exact = site,
    group = id,
    method = "rerandomize",
    options = list(p_accept = 0.2)
  )
  outer <- folds$splits[[1]]
  cal <- rsample::internal_calibration_split(
    outer,
    split_args = rsample::.get_split_args(folds)
  )
  expect_s3_class(cal, "balanced_vfold_split_cal")
  expect_s3_class(cal, "internal_calibration_split")
  expect_s3_class(cal, "balanced_vfold_split")

  fit_set <- rsample::analysis(cal)
  cal_set <- rsample::calibration(cal)
  expect_equal(nrow(fit_set) + nrow(cal_set), nrow(rsample::analysis(outer)))
  # Groups stay together, and about 1/5 of the groups go to calibration.
  expect_length(intersect(fit_set$id, cal_set$id), 0)
  n_groups <- length(unique(rsample::analysis(outer)$id))
  expect_equal(length(unique(cal_set$id)), round(n_groups / 5), tolerance = 1)

  # The inner split is balanced_initial_split() on the analysis set, with the
  # outer split's arguments.
  set.seed(30)
  cal <- rsample::internal_calibration_split(
    outer,
    split_args = rsample::.get_split_args(folds)
  )
  set.seed(30)
  direct <- balanced_initial_split(
    rsample::analysis(outer),
    prop = 1 - 1 / 5,
    balance = c(x1, x2),
    exact = site,
    group = id,
    method = "rerandomize",
    options = list(p_accept = 0.2)
  )
  expect_identical(cal$in_id, direct$in_id)
  args <- c("prop", "balance", "exact", "group", "method", "options")
  expect_identical(attributes(cal)[args], attributes(direct)[args])
})

test_that("calibration splits of an initial split use the same prop", {
  d <- sim_data()
  set.seed(4)
  split <- balanced_initial_split(d, prop = 4 / 5, balance = c(x1, cat))
  cal <- rsample::internal_calibration_split(
    split,
    split_args = rsample::.get_split_args(split)
  )
  expect_s3_class(cal, "balanced_initial_split_cal")
  expect_equal(nrow(rsample::calibration(cal)), 32)

  # tune::last_fit() passes no split arguments.
  cal <- rsample::internal_calibration_split(split, split_args = list())
  expect_equal(nrow(rsample::calibration(cal)), 32)
  expect_equal(attr(cal, "balance"), c("x1", "cat"))

  # The inner split is balanced_initial_split() on the training set, with the
  # outer split's arguments. These differ from the defaults and from the
  # fold test above, so each argument has to be passed on.
  split <- balanced_initial_split(
    d,
    prop = 2 / 3,
    balance = c(x1, cat),
    exact = site,
    group = id,
    method = "interaction",
    options = list(n_bins = 3)
  )
  set.seed(40)
  cal <- rsample::internal_calibration_split(split, split_args = list())
  set.seed(40)
  direct <- balanced_initial_split(
    rsample::training(split),
    prop = 2 / 3,
    balance = c(x1, cat),
    exact = site,
    group = id,
    method = "interaction",
    options = list(n_bins = 3)
  )
  expect_identical(cal$in_id, direct$in_id)
  args <- c("prop", "balance", "exact", "group", "method", "options")
  expect_identical(attributes(cal)[args], attributes(direct)[args])
})

test_that("twinning folds drop `strategy` for the calibration split", {
  skip_if_not_installed("twinning")
  d <- sim_data()
  set.seed(5)
  folds <- balanced_vfold_cv(
    d,
    v = 4,
    balance = c(x1, x2),
    method = "twinning",
    options = list(strategy = 2)
  )
  cal <- rsample::internal_calibration_split(
    folds$splits[[1]],
    split_args = rsample::.get_split_args(folds)
  )
  expect_equal(nrow(rsample::calibration(cal)), 150 / 4, tolerance = 1)
})

test_that("calibration falls back to an empty set when it can't split", {
  set.seed(6)
  folds <- balanced_vfold_cv(mtcars[1:4, ], v = 2, balance = mpg)
  args <- rsample::.get_split_args(folds)
  # A column that is not in the data makes the inner split fail.
  args$balance <- "not_a_column"
  expect_snapshot(
    cal <- rsample::internal_calibration_split(folds$splits[[1]], args)
  )
  expect_equal(nrow(rsample::calibration(cal)), 0)
  expect_equal(nrow(rsample::analysis(cal)), 2)
})

test_that("twinning calibration falls back to an empty set for a bad share", {
  skip_if_not_installed("twinning")
  d <- sim_data()
  set.seed(8)
  split <- balanced_initial_split(
    d,
    prop = 3 / 4,
    balance = x1,
    method = "twinning"
  )
  args <- rsample::.get_split_args(split)
  args$prop <- 0.7
  expect_snapshot(
    cal <- rsample::internal_calibration_split(split, args)
  )
  expect_equal(nrow(rsample::calibration(cal)), 0)
  expect_equal(nrow(rsample::analysis(cal)), 150)
})

test_that("tune trains a postprocessor on a balanced calibration split", {
  # tune 2.0.0 is the first version that trains postprocessors, and probably
  # 1.0.0 the first with the linear calibration.
  skip_if_not_installed("tune", "2.0.0")
  skip_if_not_installed("parsnip")
  skip_if_not_installed("workflows")
  skip_if_not_installed("tailor")
  skip_if_not_installed("probably", "1.0.0")

  # A linear calibration without smoothing fits on small calibration sets.
  # With smoothing, older probably versions skip the calibration.
  wf <- workflows::workflow(
    mpg ~ hp + wt,
    parsnip::linear_reg(),
    tailor::adjust_numeric_calibration(
      tailor::tailor(),
      method = "linear",
      smooth = FALSE
    )
  )
  # Predictions of the fitted model alone, without the postprocessor. The
  # formula uses the numeric columns as they are.
  raw_predict <- function(fitted, new_data) {
    model <- workflows::extract_fit_parsnip(fitted)
    stats::predict(model, new_data[c("hp", "wt")])$.pred
  }

  # Record the inner splits that tune gets from balancesample. All the
  # internal_calibration_split() methods make them with
  # balanced_calibration_split(). The splits themselves are tested above.
  inner <- list()
  make_inner <- balanced_calibration_split
  local_mocked_bindings(
    balanced_calibration_split = function(...) {
      split <- make_inner(...)
      inner[[length(inner) + 1]] <<- split
      split
    }
  )
  # The model in a fitted workflow was fit on the analysis set of the one
  # inner split made from `outer`, the outer analysis set.
  fit_on_inner <- function(fitted, outer) {
    from_outer <- vapply(
      inner,
      function(split) {
        # The same rows, in any order.
        rows <- rownames(outer)
        setequal(rownames(split$data), rows) &&
          identical(split$data[rows, , drop = FALSE], outer)
      },
      logical(1)
    )
    if (sum(from_outer) != 1) {
      return(FALSE)
    }
    split <- inner[[which(from_outer)]]
    refit <- stats::lm(mpg ~ hp + wt, data = rsample::analysis(split))
    engine <- workflows::extract_fit_engine(fitted)
    isTRUE(all.equal(unname(stats::coef(refit)), unname(stats::coef(engine))))
  }

  set.seed(7)
  folds <- balanced_vfold_cv(mtcars, v = 4, balance = c(mpg, hp, wt))
  ctrl <- tune::control_resamples(save_pred = TRUE, extract = function(x) x)
  res <- tune::fit_resamples(wf, resamples = folds, control = ctrl)
  expect_s3_class(res, "resample_results")
  expect_true(all(tune::collect_metrics(res)$n == 4))
  expect_length(inner, 4)
  expect_true(all(vapply(inner, inherits, logical(1), "balanced_vfold_split")))
  # In each fold, the collected predictions are the calibrated predictions
  # of the workflow that tune fitted, and the calibration changed them.
  pred <- tune::collect_predictions(res)
  extracts <- tune::collect_extracts(res)
  for (id in folds$id) {
    fitted <- extracts$.extracts[[which(extracts$id == id)]]
    fold_pred <- pred[pred$id == id, ]
    new_data <- mtcars[fold_pred$.row, ]
    expect_equal(
      unname(fold_pred$.pred),
      unname(stats::predict(fitted, new_data)$.pred)
    )
    expect_gt(max(abs(fold_pred$.pred - raw_predict(fitted, new_data))), 1e-6)
    outer <- rsample::analysis(folds$splits[[which(folds$id == id)]])
    expect_true(fit_on_inner(fitted, outer), info = id)
  }

  inner <- list()
  set.seed(8)
  split <- balanced_initial_split(mtcars, balance = c(mpg, hp, wt))
  res <- tune::last_fit(wf, split = split)
  expect_length(inner, 1)
  expect_s3_class(inner[[1]], "balanced_initial_split")
  pred <- tune::collect_predictions(res)
  expect_equal(nrow(pred), 8)
  fitted <- tune::extract_workflow(res)
  new_data <- mtcars[pred$.row, ]
  expect_equal(
    unname(pred$.pred),
    unname(stats::predict(fitted, new_data)$.pred)
  )
  expect_gt(max(abs(pred$.pred - raw_predict(fitted, new_data))), 1e-6)
  expect_true(fit_on_inner(fitted, rsample::training(split)))
})

test_that("labels() of repeated folds combine both IDs", {
  set.seed(8)
  folds <- balanced_vfold_cv(mtcars, v = 2, repeats = 2, balance = mpg)
  expect_equal(
    labels(folds),
    c("Repeat1.Fold1", "Repeat1.Fold2", "Repeat2.Fold1", "Repeat2.Fold2")
  )
  expect_s3_class(labels(folds, make_factor = TRUE), "factor")
  single <- balanced_vfold_cv(mtcars, v = 2, balance = mpg)
  expect_equal(labels(single), c("Fold1", "Fold2"))
})

test_that("calibration splits warn about their own small strata", {
  d <- mtcars
  d$s <- rep(letters[1:16], each = 2)
  set.seed(1)
  # Every stratum is in both outer sets, so the outer split is quiet.
  outer <- expect_no_warning(balanced_initial_split(d, prop = 0.51, exact = s))
  expect_warning(
    rsample::internal_calibration_split(outer, list()),
    class = "balancesample_sparse_strata"
  )
})

test_that("a small-strata warning made an error doesn't fail the split", {
  d <- mtcars
  d$s <- rep(letters[1:16], each = 2)
  set.seed(1)
  outer <- balanced_initial_split(d, prop = 0.51, exact = s)
  withr::local_options(warn = 2)
  # The error is the warning itself, not a failed calibration split.
  expect_error(
    rsample::internal_calibration_split(outer, list()),
    "^\\(converted from warning\\) 16 of 16 `exact` strata"
  )

  # The same for other warnings, such as rerandomization's.
  withr::local_options(warn = 0)
  set.seed(5)
  outer <- balanced_initial_split(
    mtcars,
    balance = c(mpg, hp, wt),
    method = "rerandomize",
    options = list(p_accept = 0.05, max_draws = 20L)
  )
  expect_warning(
    cal <- rsample::internal_calibration_split(outer, list()),
    "No draw met"
  )
  expect_equal(nrow(rsample::calibration(cal)), 6)
  withr::local_options(warn = 2)
  set.seed(5)
  outer <- balanced_initial_split(
    mtcars,
    balance = c(mpg, hp, wt),
    method = "rerandomize",
    options = list(p_accept = 0.05, max_draws = 20L)
  )
  expect_error(
    rsample::internal_calibration_split(outer, list()),
    "^\\(converted from warning\\) No draw met"
  )
})

test_that("warnings before a failed calibration split come first", {
  order <- character()
  # The inner split warns and then fails.
  local_mocked_bindings(
    balanced_initial_split = function(...) {
      warning("first")
      stop("broken")
    }
  )
  withCallingHandlers(
    cal <- balanced_calibration_split(
      mtcars,
      list(prop = 3 / 4, balance = "mpg"),
      class = c("internal_calibration_split", "rsplit")
    ),
    warning = function(cnd) {
      order <<- c(order, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(order[1], "first")
  expect_match(order[2], "Cannot create calibration split")
  expect_length(order, 2)
  expect_equal(nrow(rsample::calibration(cal)), 0)
})

test_that("calibration splits keep na = \"include\"", {
  d <- sim_data()
  d$x1[1:20] <- NA
  d$sex[21:30] <- NA

  set.seed(7)
  folds <- balanced_vfold_cv(
    d,
    v = 4,
    balance = x1,
    exact = sex,
    na = "include"
  )
  args <- rsample::.get_split_args(folds)
  expect_identical(args$na, "include")
  cal <- rsample::internal_calibration_split(folds$splits[[1]], args)
  expect_identical(attr(cal, "na"), "include")
  expect_gt(nrow(rsample::calibration(cal)), 0)

  set.seed(7)
  split <- balanced_initial_split(d, balance = x1, exact = sex, na = "include")
  # tune::last_fit() passes no split arguments.
  cal <- rsample::internal_calibration_split(split, split_args = list())
  expect_identical(attr(cal, "na"), "include")
  expect_gt(nrow(rsample::calibration(cal)), 0)

  set.seed(7)
  split3 <- balanced_initial_validation_split(
    d,
    balance = x1,
    exact = sex,
    na = "include"
  )
  cal <- rsample::internal_calibration_split(split3, split_args = list())
  expect_identical(attr(cal, "na"), "include")
  expect_gt(nrow(rsample::calibration(cal)), 0)
})

test_that("calibration splits leave out columns with no observed values", {
  d <- sim_data()
  d$z <- NA_real_
  d$z[1:5] <- 1:5
  d$s <- NA_character_
  d$s[1:5] <- "a"
  split <- rsample::make_splits(list(analysis = 6:200, assessment = 1:5), d)
  args <- list(
    balance = c("x1", "z"),
    exact = c("sex", "s"),
    na = "include",
    prop = 3 / 4
  )
  set.seed(8)
  expect_no_warning(
    cal <- balanced_calibration_split(
      rsample::analysis(split),
      args,
      class = class(split)
    )
  )
  expect_identical(attr(cal, "balance"), "x1")
  expect_identical(attr(cal, "exact"), "sex")
  expect_gt(nrow(rsample::assessment(cal)), 0)
})
