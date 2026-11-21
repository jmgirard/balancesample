# balance_check() needs columns for plain rsample objects

    Code
      balance_check(rsample::vfold_cv(d, v = 3))
    Condition
      Error in `balance_check()`:
      ! Select the columns to check.
      i For example, `balance_check(folds, age, income)`.

---

    Code
      balance_check(1)
    Condition
      Error in `balance_check()`:
      ! `x` must be an <rset>, an <rsplit>, or a <three_way_split>, not a number.

# balance_check() subsamples large sets for the energy distance

    Code
      balance_check(folds, max_n = 1)
    Condition
      Error in `balance_check()`:
      ! `max_n` must be a whole number of at least 2.

# summary() handles empty and incomplete input

    Code
      summary(bc[c("id", "statistic", "value")])
    Condition
      Error in `summary()`:
      ! `object` must have ID columns and the columns `variable`, `level`, `statistic`, and `value`.
      i Summarize the output of `balance_check()`.

---

    Code
      summary(bc[-1])
    Condition
      Error in `summary()`:
      ! `object` must have ID columns and the columns `variable`, `level`, `statistic`, and `value`.
      i Summarize the output of `balance_check()`.

---

    Code
      summary(dup)
    Condition
      Error in `summary()`:
      ! Each split in `object` must have its own ID.
      x Some IDs have more than one energy row.

---

    Code
      summary(bc[c(NA, seq_len(nrow(bc))), ])
    Condition
      Error in `summary()`:
      ! `statistic` must be "smd", "energy", or "variance_ratio" in every row of `object`.

# variance ratios show differences in spread

    Code
      balance_check(s, x, variance_ratio = NA)
    Condition
      Error in `balance_check()`:
      ! `variance_ratio` must be `TRUE` or `FALSE`.

# every split of an rset must have the columns to check

    Code
      balance_check(r, x)
    Condition
      Error in `balance_check()`:
      ! Every split of `x` must have the columns to check.
      x Split 2 has no x.

