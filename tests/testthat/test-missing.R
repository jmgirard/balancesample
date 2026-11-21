na_data <- function(n = 200) {
  d <- sim_data(n)
  withr::with_seed(7, {
    d$x1[sample(n, 30)] <- NA
    d$cat[sample(n, 20)] <- NA
    d$flag <- stats::runif(n) < 0.4
    d$flag[sample(n, 10)] <- NA
    d$sex[sample(n, 15)] <- NA
  })
  d
}

test_that("missing values are an error by default", {
  d <- na_data()
  expect_snapshot(error = TRUE, {
    balanced_vfold_cv(d, balance = c(x1, x2, cat))
    balanced_vfold_cv(d, balance = x2, na = "drop")
  })
  d$id[1:2] <- NA
  expect_snapshot(balanced_vfold_cv(d, balance = x1, group = id), error = TRUE)
})

test_that("na = \"include\" keeps the rows and spreads the missing values", {
  d <- na_data()
  for (m in available_methods()) {
    set.seed(1)
    folds <- balanced_vfold_cv(
      d,
      v = 4,
      balance = c(x1, x2, cat, flag),
      exact = sex,
      na = "include",
      method = m
    )
    expect_true(assessed_once(folds, nrow(d)), info = m)
    # The data keep their missing values.
    expect_identical(folds$splits[[1]]$data, d)
    # Missing values in `exact` are a stratum with exact counts.
    n_na <- vapply(
      folds$splits,
      function(s) sum(is.na(rsample::assessment(s)$sex)),
      integer(1)
    )
    expect_lte(diff(range(n_na)), 1)
  }
  expect_identical(attr(folds, "na"), "include")
  expect_output(print(folds), "missing values included")

  # Balancing the share missing beats random folds on it.
  share_smd <- function(folds) {
    bc <- balance_check(folds, x1, cat, flag)
    max(abs(bc$value[bc$level %in% "(missing)"]))
  }
  set.seed(2)
  random <- mean(replicate(10, share_smd(rsample::vfold_cv(d, v = 4))))
  set.seed(2)
  bal <- share_smd(balanced_vfold_cv(
    d,
    v = 4,
    balance = c(x1, cat, flag),
    na = "include"
  ))
  expect_lt(bal, random / 2)
})

test_that("na = \"include\" changes nothing without missing values", {
  d <- sim_data()
  set.seed(3)
  a <- balanced_vfold_cv(d, v = 4, balance = c(x1, cat), exact = sex)
  set.seed(3)
  b <- balanced_vfold_cv(
    d,
    v = 4,
    balance = c(x1, cat),
    exact = sex,
    na = "include"
  )
  expect_identical(
    lapply(a$splits, rsample::complement),
    lapply(b$splits, rsample::complement)
  )
})

test_that("missing values are encoded with a fill and an indicator", {
  df <- data.frame(
    x = c(1, NA, 3, 5),
    f = factor(c("a", NA, "b", "a")),
    l = c(TRUE, FALSE, NA, TRUE)
  )
  m <- build_balance_matrix(df)
  expect_identical(
    colnames(m),
    c("x", "x_(missing)", "f_a", "f_b", "f_(missing)", "l", "l_(missing)")
  )
  # The fill is the mean of the observed values, which is 0 once centered.
  expect_equal(unname(m[2, "x"]), 0)
  expect_equal(unname(m[3, "l"]), 0)
  expect_equal(m[, "x_(missing)"], m[, "f_(missing)"])

  # A real level called "NA" stays apart from missing values.
  m <- build_balance_matrix(data.frame(f = c("NA", NA, "b", "NA")))
  expect_identical(colnames(m), c("f_NA", "f_b", "f_(missing)"))
})

test_that("missing exact values are their own level in every combination", {
  d <- na_data()
  d$site2 <- d$site
  d$site2[1:12] <- NA
  set.seed(4)
  folds <- balanced_vfold_cv(
    d,
    v = 3,
    balance = x2,
    exact = c(sex, site2),
    na = "include"
  )
  cell <- interaction(
    addNA(factor(d$sex)),
    addNA(factor(d$site2)),
    drop = TRUE
  )
  for (s in folds$splits) {
    n_cell <- table(cell[rsample::complement(s)])
    share <- table(cell) / 3
    expect_true(all(n_cell >= floor(share) & n_cell <= ceiling(share)))
  }

  # Whole-number columns with missing values are still categories, and
  # infinite values are still an error.
  d$wave <- rep(1:4, length.out = nrow(d))
  d$wave[5] <- NA
  expect_no_error(
    balanced_vfold_cv(d, v = 3, balance = x2, exact = wave, na = "include")
  )
  d$wave[6] <- Inf
  expect_snapshot(
    balanced_vfold_cv(d, v = 3, balance = x2, exact = wave, na = "include"),
    error = TRUE
  )
})

test_that("na = \"include\" still rejects some missing values", {
  d <- na_data()
  d$empty <- NA_real_
  d$site[d$id == d$id[1]][1] <- NA
  expect_snapshot(error = TRUE, {
    balanced_vfold_cv(d, balance = c(x1, empty), na = "include")
    balanced_vfold_cv(d, balance = x1, exact = empty, na = "include")
    balanced_vfold_cv(d, balance = x1, exact = site, group = id, na = "include")
  })
  d$id[1:2] <- NA
  expect_snapshot(
    balanced_vfold_cv(d, balance = x1, group = id, na = "include"),
    error = TRUE
  )
})

test_that("a group with a missing exact value in every row is one stratum", {
  d <- sim_data()
  d$site[d$id %in% 1:3] <- NA
  set.seed(5)
  folds <- balanced_vfold_cv(
    d,
    v = 3,
    balance = x1,
    exact = site,
    group = id,
    na = "include"
  )
  expect_true(groups_intact(folds, d$id))
  n_groups <- vapply(
    folds$splits,
    function(s) {
      length(unique(rsample::assessment(s)$id[
        is.na(rsample::assessment(s)$site)
      ]))
    },
    integer(1)
  )
  expect_identical(sort(n_groups), c(1L, 1L, 1L))
})

test_that("the other split functions take na", {
  d <- na_data()
  set.seed(6)
  split <- balanced_initial_split(
    d,
    balance = c(x1, cat),
    exact = sex,
    na = "include"
  )
  expect_identical(attr(split, "na"), "include")
  expect_output(print(split), "missing values included")
  set.seed(6)
  split3 <- balanced_initial_validation_split(
    d,
    balance = c(x1, cat),
    exact = sex,
    na = "include"
  )
  expect_identical(attr(split3, "val_att")$na, "include")
  expect_output(print(split3), "missing values included")
  expect_snapshot(
    balanced_initial_split(d, balance = c(x1, cat)),
    error = TRUE
  )
})
