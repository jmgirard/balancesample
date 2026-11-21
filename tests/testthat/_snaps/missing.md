# missing values are an error by default

    Code
      balanced_vfold_cv(d, balance = c(x1, x2, cat))
    Condition
      Error in `balanced_vfold_cv()`:
      ! Selected columns must not have missing values.
      x Missing values in x1 (30 rows) and cat (20 rows).
      i Use `na = "include"` to balance on which values are missing, or drop those rows. See `vignette("balancesample")`.
    Code
      balanced_vfold_cv(d, balance = x2, na = "drop")
    Condition
      Error in `balanced_vfold_cv()`:
      ! `na` must be one of "error" or "include", not "drop".

---

    Code
      balanced_vfold_cv(d, balance = x1, group = id)
    Condition
      Error in `balanced_vfold_cv()`:
      ! Selected columns must not have missing values.
      x Missing values in x1 (30 rows) and id (2 rows).
      i Use `na = "include"` to balance on which values are missing, or drop those rows. See `vignette("balancesample")`.
      i `group` can't have missing values. Give the rows with a missing ID their own IDs, or one shared ID, depending on whether they can belong together.

# missing exact values are their own level in every combination

    Code
      balanced_vfold_cv(d, v = 3, balance = x2, exact = wave, na = "include")
    Condition
      Error in `balanced_vfold_cv()`:
      ! wave has infinite values.

# na = "include" still rejects some missing values

    Code
      balanced_vfold_cv(d, balance = c(x1, empty), na = "include")
    Condition
      Error in `balanced_vfold_cv()`:
      ! `balance` and `exact` columns must have at least one observed value.
      x empty is missing in every row.
    Code
      balanced_vfold_cv(d, balance = x1, exact = empty, na = "include")
    Condition
      Error in `balanced_vfold_cv()`:
      ! `balance` and `exact` columns must have at least one observed value.
      x empty is missing in every row.
    Code
      balanced_vfold_cv(d, balance = x1, exact = site, group = id, na = "include")
    Condition
      Error in `balanced_vfold_cv()`:
      ! `exact` columns must be constant within each group.
      x 1 group contains more than one `exact` level.
      i A missing value counts as its own level. Fill it in where you know the group's value.

---

    Code
      balanced_vfold_cv(d, balance = x1, group = id, na = "include")
    Condition
      Error in `balanced_vfold_cv()`:
      ! `group` must not have missing values.
      x id is missing in 2 rows.
      i Give the rows with a missing ID their own IDs, or one shared ID, depending on whether they can belong together.

# the other split functions take na

    Code
      balanced_initial_split(d, balance = c(x1, cat))
    Condition
      Error in `balanced_initial_split()`:
      ! Selected columns must not have missing values.
      x Missing values in x1 (30 rows) and cat (20 rows).
      i Use `na = "include"` to balance on which values are missing, or drop those rows. See `vignette("balancesample")`.

