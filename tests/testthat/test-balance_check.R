test_that("balance_check() returns one SMD per variable and level", {
  d <- sim_data()
  set.seed(1)
  folds <- balanced_vfold_cv(
    d,
    v = 4,
    balance = c(x1, cat),
    exact = sex,
    method = "rerandomize"
  )
  bc <- balance_check(folds)
  expect_s3_class(bc, "balance_check")
  expect_named(bc, c("id", "variable", "level", "statistic", "value"))
  # x1 + 3 levels of cat + 2 levels of sex + energy = 7 rows per fold
  expect_equal(nrow(bc), 4 * 7)
  expect_equal(sum(bc$statistic == "energy"), 4)
  expect_true(all(bc$value[bc$statistic == "energy"] >= 0))
})

test_that("balance_check() matches a direct SMD", {
  d <- sim_data()
  set.seed(2)
  folds <- rsample::vfold_cv(d, v = 3)
  bc <- balance_check(folds, x1)
  s <- folds$splits[[1]]
  a <- d$x1[rsample::complement(s)]
  b <- d$x1[s$in_id]
  expected <- (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)
  expect_equal(bc$value[bc$id == "Fold1" & bc$statistic == "smd"], expected)
})

test_that("energy distance is near zero for identical sets and grows with a shift", {
  set.seed(3)
  a <- matrix(rnorm(200), ncol = 2)
  expect_equal(energy_distance(a, a), 0)
  expect_gt(energy_distance(a, a + 2), energy_distance(a, a + 0.1))
})

test_that("balance_check() works on repeats and on splits", {
  d <- sim_data()
  set.seed(4)
  folds <- balanced_vfold_cv(
    d,
    v = 3,
    repeats = 2,
    balance = x1,
    method = "rerandomize"
  )
  bc <- balance_check(folds)
  expect_true(all(c("id", "id2") %in% names(bc)))

  s <- balanced_initial_split(d, balance = c(x1, x2), method = "rerandomize")
  bc <- balance_check(s)
  expect_equal(unique(bc$id), "Split")
  expect_output(print(bc), "Balance check: 1 split")
  withr::local_options(width = 40)
  header <- utils::capture.output(print(bc))[1:3]
  expect_match(header[1], "^# Balance check: 1 split,$")
  expect_match(header[2], "^# max \\|SMD\\| = [0-9.e-]+,$")
  expect_match(header[3], "^# mean energy distance = [0-9.e-]+$")

  withr::local_options(width = 80)
  full <- utils::capture.output(print(bc))[1]
  withr::local_options(width = nchar(full))
  expect_equal(utils::capture.output(print(bc))[1], full)
  withr::local_options(width = nchar(full) - 1)
  expect_match(utils::capture.output(print(bc))[1], ",$")
})

test_that("balance_check() drops missing values", {
  d <- data.frame(z = c(1, 2, 3, NA, 5, 6), w = c("a", NA, "b", "a", "b", "a"))
  s <- rsample::make_splits(list(analysis = 1:3, assessment = 4:6), d)
  bc <- balance_check(s, z, w)
  expect_equal(bc$value[bc$variable %in% "z"], (5.5 - 2) / sqrt((0.5 + 1) / 2))
  expect_false(anyNA(bc$value))
})

test_that("balance_check() reports complete separation as Inf", {
  d <- data.frame(g = rep(c("A", "B"), each = 4))
  s <- rsample::make_splits(list(analysis = 1:4, assessment = 5:8), d)
  bc <- balance_check(s, g)
  expect_equal(bc$value[bc$statistic == "smd"], c(-Inf, Inf))
  expect_output(print(bc), "max \\|SMD\\| = Inf")

  d <- data.frame(x = c(0, 0, 1e-13, 1e-13))
  s <- rsample::make_splits(list(analysis = 1:2, assessment = 3:4), d)
  expect_equal(balance_check(s, x)$value[1], Inf)
})

test_that("balance_check() needs columns for plain rsample objects", {
  d <- sim_data()
  expect_snapshot(balance_check(rsample::vfold_cv(d, v = 3)), error = TRUE)
  expect_snapshot(balance_check(1), error = TRUE)
})

test_that("autoplot() returns ggplot objects", {
  skip_if_not_installed("ggplot2")
  d <- sim_data()
  set.seed(5)
  folds <- balanced_vfold_cv(
    d,
    v = 3,
    repeats = 2,
    balance = c(x1, cat),
    method = "rerandomize"
  )
  bc <- balance_check(folds)
  expect_s3_class(ggplot2::autoplot(bc), "ggplot")
  expect_s3_class(ggplot2::autoplot(bc, type = "energy"), "ggplot")
  # The record of ID columns doesn't reach the plotted data.
  expect_null(attr(ggplot2::autoplot(bc)$data, "id_cols"))

  s <- balanced_initial_split(d, balance = x1, method = "rerandomize")
  p <- ggplot2::autoplot(balance_check(s), type = "energy")
  expect_s3_class(p, "ggplot")
})

test_that("balance_check() handles logical and missing columns", {
  d <- data.frame(
    flag = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE),
    z = c(NA, NA, NA, 1, 2, 3)
  )
  s <- rsample::make_splits(list(analysis = 1:3, assessment = 4:6), d)
  bc <- balance_check(s, flag, z)
  expect_equal(bc$level[bc$variable %in% "flag"], c("FALSE", "TRUE"))
  # z is missing in the whole analysis set, so no row there is complete.
  expect_true(is.na(bc$value[bc$variable %in% "z"]))
  expect_true(is.na(bc$value[bc$statistic == "energy"]))
})

test_that("balance_check() subsamples large sets for the energy distance", {
  a <- matrix(stats::rnorm(60), ncol = 2)
  b <- matrix(stats::rnorm(80, mean = 1), ncol = 2)
  set.seed(6)
  sub <- energy_distance(a, b, max_n = 10)
  set.seed(6)
  a_rows <- sample.int(nrow(a), 10)
  b_rows <- sample.int(nrow(b), 10)
  expect_equal(sub, energy_distance(a[a_rows, ], b[b_rows, ]))

  d <- sim_data()
  set.seed(6)
  folds <- balanced_vfold_cv(d, v = 2, balance = c(x1, x2))
  # A smaller `max_n` changes the energy distance through the public API.
  energy <- function(bc) bc$value[bc$statistic == "energy"]
  set.seed(7)
  small <- energy(balance_check(folds, max_n = 5))
  full <- energy(balance_check(folds))
  expect_false(isTRUE(all.equal(small, full)))
  expect_error(balance_check(folds, max_n = 1), "at least 2")
  expect_snapshot(balance_check(folds, max_n = 1), error = TRUE)
})

test_that("balance_check() needs at least one split", {
  empty <- rsample::manual_rset(list(), character())
  expect_error(balance_check(empty, mpg), "no splits", class = "rlang_error")
})

test_that("summary() gives one row per split", {
  d <- sim_data()
  set.seed(8)
  folds <- balanced_vfold_cv(
    d,
    v = 3,
    repeats = 2,
    balance = c(x1, cat),
    exact = sex,
    method = "rerandomize"
  )
  bc <- balance_check(folds)
  s <- summary(bc)
  expect_false(inherits(s, "balance_check"))
  expect_named(
    s,
    c("id", "id2", "max_abs_smd", "variable", "level", "energy")
  )
  expect_equal(s$id, folds$id)
  expect_equal(s$id2, folds$id2)
  expect_equal(s$energy, bc$value[bc$statistic == "energy"])
  for (i in seq_len(nrow(s))) {
    rows <- bc[bc$id == s$id[i] & bc$id2 == s$id2[i], ]
    smd <- rows[rows$statistic == "smd", ]
    worst <- which.max(abs(smd$value))
    expect_equal(s$max_abs_smd[i], abs(smd$value[worst]))
    expect_equal(s$variable[i], smd$variable[worst])
    expect_equal(s$level[i], smd$level[worst])
  }
  expect_error(summary(bc, 1), class = "rlib_error_dots_nonempty")
})

test_that("summary() handles splits, three-way splits, and missing SMDs", {
  d <- sim_data()
  set.seed(9)
  s3 <- balanced_initial_validation_split(
    d,
    prop = c(0.6, 0.2),
    balance = c(x1, x2),
    method = "rerandomize"
  )
  bc3 <- balance_check(s3)
  out3 <- summary(bc3)
  expect_equal(out3$id, c("Validation", "Testing"))
  expect_equal(out3$energy, bc3$value[bc3$statistic == "energy"])
  for (id in out3$id) {
    smd <- bc3$value[bc3$id == id & bc3$statistic == "smd"]
    expect_equal(out3$max_abs_smd[out3$id == id], max(abs(smd)))
  }

  d <- data.frame(
    z = c(NA, NA, NA, 1, 2, 3),
    g = c("a", "b", "a", "b", "a", "b")
  )
  s <- rsample::make_splits(list(analysis = 1:3, assessment = 4:6), d)
  out <- summary(balance_check(s, z))
  expect_equal(out$id, "Split")
  expect_true(is.na(out$max_abs_smd))
  expect_true(is.na(out$variable))
  expect_true(is.na(out$energy))

  out <- summary(balance_check(s, z, g))
  expect_equal(out$variable, "g")
  # The two levels of g tie, and the first one is kept.
  expect_equal(out$level, "a")
  expect_equal(out$max_abs_smd, abs(balance_check(s, g)$value[1]))
  d$g <- factor(d$g, levels = c("b", "a"))
  s <- rsample::make_splits(list(analysis = 1:3, assessment = 4:6), d)
  expect_equal(summary(balance_check(s, g))$level, "b")

  d <- data.frame(g = rep(c("A", "B"), each = 4))
  s <- rsample::make_splits(list(analysis = 1:4, assessment = 5:8), d)
  out <- summary(balance_check(s, g))
  expect_equal(out$max_abs_smd, Inf)
  expect_equal(out$level, "A")
})

test_that("summary() keeps splits apart however the rows are arranged", {
  d <- sim_data()
  set.seed(10)
  folds <- rsample::vfold_cv(d, v = 3)
  bc <- balance_check(folds, x1, cat)
  shuffled <- bc[sample(nrow(bc)), ]
  expect_equal(
    summary(shuffled)[order(summary(shuffled)$id), ],
    summary(bc)
  )

  # IDs that would collide if pasted, and `NA` next to "NA".
  folds <- rsample::vfold_cv(d, v = 2, repeats = 2)
  bc <- balance_check(folds, x1, cat)
  key <- paste(bc$id, bc$id2)
  fold <- match(key, unique(key))
  bc$id <- c("a\rb", "a", NA, "NA")[fold]
  bc$id2 <- c("c", "b\rc", "x", "x")[fold]
  out <- summary(bc)
  expect_equal(out$id, c("a\rb", "a", NA, "NA"))
  expect_equal(out$energy, bc$value[bc$statistic == "energy"])
})

test_that("summary() handles empty and incomplete input", {
  d <- sim_data()
  set.seed(11)
  bc <- balance_check(rsample::vfold_cv(d, v = 3), x1)
  empty <- summary(bc[0, ])
  expect_equal(nrow(empty), 0)
  expect_named(empty, names(summary(bc)))
  expect_type(empty$max_abs_smd, "double")
  expect_type(empty$variable, "character")
  expect_type(empty$energy, "double")

  expect_snapshot(summary(bc[c("id", "statistic", "value")]), error = TRUE)
  expect_snapshot(summary(bc[-1]), error = TRUE)

  dup <- bc
  dup$id <- "same"
  expect_snapshot(summary(dup), error = TRUE)

  # Only the ID columns balance_check() made are IDs, even if an added
  # column's name starts with "id".
  noted <- bc
  noted$note <- seq_len(nrow(bc))
  noted$identifier <- seq_len(nrow(bc))
  expect_equal(summary(noted), summary(bc))

  # Without the record, columns whose names start with "id" are IDs.
  plain <- bc
  attr(plain, "id_cols") <- NULL
  plain$note <- seq_len(nrow(bc))
  expect_equal(summary(plain), summary(bc))
  plain$identifier <- seq_len(nrow(bc))
  expect_equal(nrow(summary(plain)), nrow(bc))

  # A renamed ID column is still an ID.
  set.seed(12)
  rep_bc <- balance_check(rsample::vfold_cv(d, v = 3, repeats = 2), x1)
  renamed <- rep_bc
  names(renamed)[names(renamed) == "id2"] <- "id_fold"
  out <- summary(renamed[renamed$statistic == "smd", ])
  expect_equal(nrow(out), 6)
  expect_equal(out$id_fold, summary(rep_bc)$id2)

  # print() returns its input, and summary() drops the record.
  expect_output(printed <- print(bc))
  expect_identical(printed, bc)
  expect_null(attr(summary(bc), "id_cols"))

  # Rows with a missing or unknown statistic.
  expect_snapshot(summary(bc[c(NA, seq_len(nrow(bc))), ]), error = TRUE)
})
