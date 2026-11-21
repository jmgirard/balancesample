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

test_that("balance_check() drops missing values and reports the share missing", {
  d <- data.frame(z = c(1, 2, 3, NA, 5, 6), w = c("a", NA, "b", "a", "b", "a"))
  s <- rsample::make_splits(list(analysis = 1:3, assessment = 4:6), d)
  bc <- balance_check(s, z, w)
  smd <- bc$statistic == "smd"
  expect_equal(
    bc$value[smd & bc$variable %in% "z" & is.na(bc$level)],
    (5.5 - 2) / sqrt((0.5 + 1) / 2)
  )
  expect_equal(bc$level[bc$variable %in% "w"], c("a", "b", "(missing)"))
  # z is missing in 1 of 3 assessment rows and in no analysis row; w the
  # other way round.
  missing <- bc$value[bc$level %in% "(missing)"]
  expect_equal(missing, c(1, -1))
  expect_false(anyNA(bc$value))
})

test_that("the energy distance uses all rows when values are missing", {
  set.seed(1)
  n <- 200
  d <- data.frame(x = rnorm(n), y = rnorm(n))
  # Every row has a missing value, so no row is complete.
  d$x[seq(1, n, 2)] <- NA
  d$y[seq(2, n, 2)] <- NA
  split <- rsample::initial_split(d)
  en <- balance_check(split, x, y)
  en <- en$value[en$statistic == "energy"]
  expect_true(is.finite(en) && en > 0)

  # The same encoding as the split functions: filled values and indicators.
  xs <- build_balance_matrix(d)
  expect_equal(
    en,
    energy_distance(
      xs[rsample::complement(split), ],
      xs[split$in_id, ]
    )
  )

  # With no missing values, nothing changes.
  d <- data.frame(x = rnorm(n), g = sample(letters[1:3], n, TRUE))
  split <- rsample::initial_split(d)
  bc <- balance_check(split, x, g)
  xs <- build_balance_matrix(d)
  expect_equal(
    bc$value[bc$statistic == "energy"],
    energy_distance(xs[rsample::complement(split), ], xs[split$in_id, ])
  )
  expect_false(any(bc$level %in% "(missing)"))
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

test_that("SMDs don't depend on the units or the origin", {
  x <- c(-2, -1, 1, 2, -1.5, 0.5)
  smd <- function(v) {
    s <- rsample::make_splits(
      list(analysis = 1:3, assessment = 4:6),
      data.frame(v = v)
    )
    balance_check(s, v)$value[1]
  }
  expect_equal(smd(x * 1e-14), smd(x))
  expect_equal(smd(x * 1e14), smd(x))
  expect_equal(smd(x + 1e12), smd(x))
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
  # z is missing in the whole analysis set, so it has no SMD of its values,
  # and its share missing differs completely.
  z <- bc$variable %in% "z"
  expect_true(is.na(bc$value[z & is.na(bc$level)]))
  expect_equal(bc$value[z & bc$level %in% "(missing)"], -Inf)
  expect_true(is.finite(bc$value[bc$statistic == "energy"]))
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
  # z has no observed value in the analysis set, so its only SMD that
  # isn't missing is the one of its share missing.
  bc <- balance_check(s, z)
  out <- summary(bc[!bc$level %in% "(missing)", ])
  expect_equal(out$id, "Split")
  expect_true(is.na(out$max_abs_smd))
  expect_true(is.na(out$variable))
  out <- summary(bc)
  expect_equal(out$variable, "z")
  expect_equal(out$level, "(missing)")
  expect_equal(out$max_abs_smd, Inf)
  bc$value[bc$statistic == "energy"] <- NA
  expect_true(is.na(summary(bc)$energy))

  out <- summary(balance_check(s, g))
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

test_that("whole-number exact columns are checked as categories", {
  d <- data.frame(code = rep(1:3, 20), x = seq_len(60))
  set.seed(1)
  folds <- balanced_vfold_cv(d, v = 3, balance = x, exact = code)
  bc <- balance_check(folds)
  expect_equal(bc$level[bc$variable %in% "code"], rep(c("1", "2", "3"), 3))

  # A set with only code 2 differs from one with codes 1 and 3, although
  # the means are equal.
  d <- data.frame(code = c(2, 2, 2, 2, 1, 3, 1, 3), y = 1:8)
  plain <- rsample::make_splits(list(analysis = 5:8, assessment = 1:4), d)
  # Named columns from a plain split stay numeric.
  expect_equal(balance_check(plain, code)$value[1], 0)
  stored <- plain
  attr(stored, "exact") <- "code"
  bc <- balance_check(stored)
  expect_equal(bc$value[bc$variable %in% "code" & bc$level == "2"], Inf)
  # Naming a stored exact column keeps it categorical, with no variance
  # ratio, and the energy distance matches that of a factor.
  bc <- balance_check(stored, code, variance_ratio = TRUE)
  expect_equal(bc$level[bc$statistic == "smd"], c("1", "2", "3"))
  expect_false(any(bc$statistic == "variance_ratio"))
  as_factor <- rsample::make_splits(
    list(analysis = 5:8, assessment = 1:4),
    transform(d, code = factor(code))
  )
  expect_equal(
    bc$value[bc$statistic == "energy"],
    balance_check(as_factor, code)$value[4]
  )
})

test_that("the SMD plot of a three-way split tells the sets apart", {
  skip_if_not_installed("ggplot2")
  d <- sim_data()
  set.seed(3)
  s <- balanced_initial_validation_split(d, balance = c(x1, x2))
  built <- ggplot2::ggplot_build(ggplot2::autoplot(balance_check(s)))
  points <- built$data[[length(built$data)]]
  expect_equal(length(unique(points$shape)), 2)
  folds <- balanced_vfold_cv(d, v = 3, balance = c(x1, x2))
  built <- ggplot2::ggplot_build(ggplot2::autoplot(balance_check(folds)))
  expect_equal(length(unique(built$data[[length(built$data)]]$shape)), 1)
})

test_that("variance ratios show differences in spread", {
  d <- data.frame(
    x = c(-1, 1, -1, 1, -5, 5, -5, 5),
    g = rep(c("a", "b"), 4),
    k = c(1, 1, 1, 1, 2, 2, 2, 2)
  )
  s <- rsample::make_splits(list(analysis = 1:4, assessment = 5:8), d)
  bc <- balance_check(s, x, g, k, variance_ratio = TRUE)
  vr <- bc[bc$statistic == "variance_ratio", ]
  # Only numeric variables get one.
  expect_equal(vr$variable, c("x", "k"))
  expect_equal(bc$value[bc$variable %in% "x" & bc$statistic == "smd"], 0)
  expect_equal(vr$value[1], 25)
  # Two zero variances give 1.
  expect_equal(vr$value[2], 1)
  expect_equal(summary(bc)$variance_ratio, 25)
  expect_false("variance_ratio" %in% names(summary(balance_check(s, x))))

  # A zero variance on one side gives 0 or Inf, and a set of one gives NA.
  d <- data.frame(x = c(1, 2, 3, 4, 4, 4))
  s <- rsample::make_splits(list(analysis = 1:3, assessment = 4:6), d)
  bc <- balance_check(s, x, variance_ratio = TRUE)
  expect_equal(bc$value[bc$statistic == "variance_ratio"], 0)
  s <- rsample::make_splits(list(analysis = 4:6, assessment = 1:3), d)
  bc <- balance_check(s, x, variance_ratio = TRUE)
  expect_equal(bc$value[bc$statistic == "variance_ratio"], Inf)
  s <- rsample::make_splits(list(analysis = 1:5, assessment = 6L), d)
  bc <- balance_check(s, x, variance_ratio = TRUE)
  expect_true(is.na(bc$value[bc$statistic == "variance_ratio"]))
  expect_true(is.na(summary(bc)$variance_ratio))

  # The ratio farthest from 1 on the log scale is summarized.
  d <- data.frame(
    x = c(-2, 2, -2, 2, -1, 1, -1, 1),
    y = c(-1, 1, -1, 1, -1.5, 1.5, -1.5, 1.5)
  )
  s <- rsample::make_splits(list(analysis = 1:4, assessment = 5:8), d)
  expect_equal(
    summary(balance_check(s, x, y, variance_ratio = TRUE))$variance_ratio,
    0.25
  )

  expect_snapshot(balance_check(s, x, variance_ratio = NA), error = TRUE)
})

test_that("variance ratios work for every kind of object", {
  d <- sim_data()
  set.seed(4)
  folds <- balanced_vfold_cv(d, v = 3, balance = c(x1, cat))
  bc <- balance_check(folds, variance_ratio = TRUE)
  expect_equal(sum(bc$statistic == "variance_ratio"), 3)
  s3 <- balanced_initial_validation_split(d, balance = x1)
  bc <- balance_check(s3, variance_ratio = TRUE)
  expect_equal(
    bc$id[bc$statistic == "variance_ratio"],
    c("Validation", "Testing")
  )
  expect_output(print(bc), "Balance check: 2 splits")
  skip_if_not_installed("ggplot2")
  expect_s3_class(ggplot2::autoplot(bc), "ggplot")
})

test_that("each split of an rset is checked on its own data", {
  a <- rsample::make_splits(
    list(analysis = 1:4, assessment = 5:8),
    data.frame(x = rep(0:1, each = 4))
  )
  b <- rsample::make_splits(
    list(analysis = 1:4, assessment = 5:8),
    data.frame(x = rep(0:1, 4))
  )
  r <- rsample::manual_rset(list(a, b), c("First", "Second"))
  s <- summary(balance_check(r, x))
  expect_equal(s$max_abs_smd, c(Inf, 0))
  expect_equal(s$energy[2], 0)
})

test_that("small differences far from the origin are measured", {
  t0 <- as.POSIXct("2026-01-01", tz = "UTC")
  # Microseconds on a date-time in seconds.
  d <- data.frame(t = t0 + c(0, 10, 0, 10, 10, 30, 10, 30) / 1e6)
  s <- rsample::make_splits(list(analysis = 1:4, assessment = 5:8), d)
  bc <- balance_check(s, t, variance_ratio = TRUE)
  d0 <- data.frame(t = c(0, 10, 0, 10, 10, 30, 10, 30))
  s0 <- rsample::make_splits(list(analysis = 1:4, assessment = 5:8), d0)
  expect_equal(bc$value, balance_check(s0, t, variance_ratio = TRUE)$value)
  expect_equal(ncol(build_balance_matrix(d)), 1)
})

test_that("a gap between the sets doesn't hide the spread within them", {
  # Spreads of 1e-7 seconds at the epoch, against a set in 2026.
  d <- data.frame(
    t = as.POSIXct(
      c(0, 1e-7, 2e-7, 0, rep(1767225600, 4)),
      origin = "1970-01-01",
      tz = "UTC"
    )
  )
  x <- as.numeric(d$t)
  # The 2026 set has no spread.
  smd <- (x[5] - mean(x[1:4])) / sqrt(var(x[1:4]) / 2)
  s <- rsample::make_splits(list(analysis = 1:4, assessment = 5:8), d)
  expect_equal(balance_check(s, t)$value[1], smd)
  # Swapping the sets flips the sign but keeps the size.
  r <- rsample::make_splits(list(analysis = 5:8, assessment = 1:4), d)
  expect_equal(balance_check(r, t)$value[1], -smd)

  # Zero spread on one side gives a variance ratio of 0.
  bc <- balance_check(s, t, variance_ratio = TRUE)
  expect_equal(bc$value[bc$statistic == "variance_ratio"], 0)
})

test_that("every split of an rset must have the columns to check", {
  a <- rsample::make_splits(
    list(analysis = 1:4, assessment = 5:8),
    data.frame(x = 1:8, y = 1:8)
  )
  b <- rsample::make_splits(
    list(analysis = 1:4, assessment = 5:8),
    data.frame(y = 1:8)
  )
  r <- rsample::manual_rset(list(a, b), c("A", "B"))
  expect_snapshot(balance_check(r, x), error = TRUE)
})

test_that("logical exact columns are checked as categories", {
  set.seed(3)
  d <- data.frame(
    flag = rep(c(FALSE, TRUE, NA), c(11, 12, 10)),
    x = rnorm(33)
  )
  split <- balanced_initial_split(d, balance = x, exact = flag, na = "include")
  bc <- balance_check(split)
  expect_equal(
    bc$level[bc$variable %in% "flag"],
    c("FALSE", "TRUE", "(missing)")
  )
  as_factor <- rsample::make_splits(
    list(analysis = split$in_id, assessment = rsample::complement(split)),
    transform(d, flag = factor(flag))
  )
  energy <- function(bc) bc$value[bc$statistic == "energy"]
  expect_equal(energy(bc), energy(balance_check(as_factor, x, flag)))
})
