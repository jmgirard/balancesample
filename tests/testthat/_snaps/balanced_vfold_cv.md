# exact columns must be constant within groups

    Code
      balanced_vfold_cv(d, v = 4, balance = x1, exact = sex, group = id, method = "rerandomize")
    Condition
      Error in `balanced_vfold_cv()`:
      ! `exact` columns must be constant within each group.
      x 35 groups contain more than one `exact` level.

# method options are checked

    Code
      balanced_vfold_cv(d, balance = x1, method = "rerandomize", options = list(foo = 1))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Unknown option for method "rerandomize": `foo`.
      i Valid options: `p_accept` and `max_draws`.
    Code
      balanced_vfold_cv(d, balance = x1, method = "nope")
    Condition
      Error in `balanced_vfold_cv()`:
      ! `method` must be one of "anticlust", "interaction", "rerandomize", or "twinning", not "nope".
    Code
      balanced_vfold_cv(d, balance = x1, method = "interaction", options = list(1))
    Condition
      Error in `balanced_vfold_cv()`:
      ! All elements of `options` must have unique names.

# anticlust options work and are checked

    Code
      balanced_vfold_cv(d, balance = x1, options = list(moments = 5))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "anticlust" failed.
      Caused by error:
      ! `moments` must be 1, 2, 3, or 4.
    Code
      balanced_vfold_cv(d, balance = x1, options = list(k_neighbours = 0.5))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "anticlust" failed.
      Caused by error:
      ! `k_neighbours` must be a whole number of at least 1, or `Inf`.
    Code
      balanced_vfold_cv(d, balance = x1, options = list(objective = "diversity"))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Unknown option for method "anticlust": `objective`.
      i Valid options: `moments`, `covariances`, `k_neighbours`, `local_maximum`, and `repetitions`.

# twinning strategy 2 needs a power of two

    Code
      balanced_vfold_cv(d, v = 5, balance = x1, method = "twinning", options = list(
        strategy = 2))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "twinning" failed.
      Caused by error:
      ! Twinning strategy 2 needs `v` to be a power of 2.
      x `v` is 5.
      i Use `v` = 2, 4, 8, or 16, or use strategy 1 or 3.

# bad inputs give clear errors

    Code
      balanced_vfold_cv(d)
    Condition
      Error in `balanced_vfold_cv()`:
      ! Select at least one column with `balance` or `exact`.
      i For example, `balance = c(age, income)`.
    Code
      balanced_vfold_cv(d, v = 1, balance = x1)
    Condition
      Error in `balanced_vfold_cv()`:
      ! `v` must be a whole number of at least 2.
    Code
      balanced_vfold_cv(d[1:3, ], v = 4, balance = x1)
    Condition
      Error in `balanced_vfold_cv()`:
      ! `v` (4) must not be more than the number of rows (3).
    Code
      balanced_vfold_cv(d, balance = x1, exact = x2)
    Condition
      Error in `balanced_vfold_cv()`:
      ! `exact` columns must be categorical.
      x x2 is a double vector.
      i Convert it to a factor, or move it to `balance`.
    Code
      balanced_vfold_cv(d, balance = c(x1, id), group = id)
    Condition
      Error in `balanced_vfold_cv()`:
      ! The `group` column id can't also be in `balance` or `exact`.
    Code
      balanced_vfold_cv(d, balance = x1, group = c(id, sex))
    Condition
      Error in `balanced_vfold_cv()`:
      ! `group` must select one column, not 2.
    Code
      balanced_vfold_cv(d, v = 50, balance = x1, group = id)
    Condition
      Error in `balanced_vfold_cv()`:
      ! `v` (50) must not be more than the number of groups (38).

---

    Code
      balanced_vfold_cv(d, balance = x1)
    Condition
      Error in `balanced_vfold_cv()`:
      ! Selected columns must not have missing values.
      x Missing values in x1.

---

    Code
      balanced_vfold_cv(m, 3, balance = x)
    Condition
      Error in `balanced_vfold_cv()`:
      ! Selected columns must be plain vectors.
      x x is a matrix, data frame, or list column.

# strata and other extra arguments are rejected

    Code
      balanced_vfold_cv(d, balance = x1, strata = cat)
    Condition
      Error in `balanced_vfold_cv()`:
      ! balancesample uses `exact`, not `strata`.
      i Put categorical columns in `exact`, for example `exact = class`.
      i Put numeric columns in `balance`.
    Code
      balanced_initial_split(d, balance = x1, strata = cat)
    Condition
      Error in `balanced_initial_split()`:
      ! balancesample uses `exact`, not `strata`.
      i Put categorical columns in `exact`, for example `exact = class`.
      i Put numeric columns in `balance`.
    Code
      balanced_vfold_cv(d, balance = x1, folds = 5)
    Condition
      Error in `balanced_vfold_cv()`:
      ! `...` must be empty.
      x Problematic argument:
      * folds = 5

# bad data give clear errors

    Code
      balanced_vfold_cv(1:10, balance = x)
    Condition
      Error in `balanced_vfold_cv()`:
      ! `data` must be a data frame, not an integer vector.
    Code
      balanced_vfold_cv(data.frame(x = 1), balance = x)
    Condition
      Error in `balanced_vfold_cv()`:
      ! `data` must have at least two rows.
    Code
      balanced_vfold_cv(data.frame(x = c(1, Inf, 3, 4)), v = 2, balance = x)
    Condition
      Error in `balanced_vfold_cv()`:
      ! x has infinite values.
    Code
      balanced_vfold_cv(data.frame(x = as.Date("2026-01-01") + 0:3), v = 2, balance = x)
    Condition
      Error in `balanced_vfold_cv()`:
      ! `balance` columns must be numeric, logical, character, or factor.
      x x is a <Date> object.

# more anticlust options are checked

    Code
      balanced_vfold_cv(d, balance = x1, options = list(covariances = "yes"))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "anticlust" failed.
      Caused by error:
      ! `covariances` must be `TRUE` or `FALSE`.
    Code
      balanced_vfold_cv(d, balance = x1, options = list(local_maximum = NA))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "anticlust" failed.
      Caused by error:
      ! `local_maximum` must be `TRUE` or `FALSE`.
    Code
      balanced_vfold_cv(d, balance = x1, options = list(repetitions = 0))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "anticlust" failed.
      Caused by error:
      ! `repetitions` must be a whole number of at least 1.

# rerandomize and interaction options are checked

    Code
      balanced_vfold_cv(d, balance = x1, method = "rerandomize", options = list(
        p_accept = 0))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "rerandomize" failed.
      Caused by error:
      ! `p_accept` must be a number above 0 and at most 1.
    Code
      balanced_vfold_cv(d, balance = x1, method = "rerandomize", options = list(
        max_draws = 0))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "rerandomize" failed.
      Caused by error:
      ! `max_draws` must be a whole number of at least 1.
    Code
      balanced_vfold_cv(d, balance = x1, method = "interaction", options = list(
        n_bins = 1))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "interaction" failed.
      Caused by error:
      ! `n_bins` must be a whole number of at least 2.
    Code
      balanced_vfold_cv(d, balance = x1, method = "interaction", options = list(pool = 1))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Method "interaction" failed.
      Caused by error:
      ! `pool` must be a number from 0 up to (not including) 1.

