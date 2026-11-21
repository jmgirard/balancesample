# balancesample 0.1.0

* First CRAN release.

* `balanced_vfold_cv()`, `balanced_initial_split()`, and
  `balanced_initial_validation_split()` make v-fold cross-validation folds,
  training/testing splits, and training/validation/testing splits that are
  balanced on several variables at once. The three sets of a three-way split
  are balanced in one pass. Columns in `exact` are stratified exactly: each
  combination of their levels is a stratum, and its count in each set is its
  share of that set's actual size, rounded up or down. `group` keeps the rows
  of a group in the same set; then the proportions and the exact counts are
  numbers of groups, and `exact` columns must be constant within each group.
  If more than half of the `exact` strata are missing from at least one
  set, you get a warning. `balance` takes numeric, logical, character,
  factor, `Date`, `POSIXct`, and `difftime` columns. Printed splits and
  folds show the method and the columns.

* By default, the selected columns must not have missing values. With
  `na = "include"`, rows with missing values in `balance` and `exact`
  columns stay in, and the split also balances which values are missing. A
  numeric column gets the mean of its observed values plus a 0/1 column of
  which values are missing, and in a categorical or `exact` column a
  missing value is one more level. The split's data keep their missing
  values. `group` must not have missing values either way. Inner
  calibration splits made by tune use the same `na`.

* Four balancing methods: `"anticlust"` (k-plus anticlustering, the
  default), `"twinning"` (data twinning; needs the twinning package),
  `"rerandomize"` (rerandomization on the Mahalanobis imbalance), and
  `"interaction"` (random stratification on crossed quantile bins, as a
  baseline). `?balanced_vfold_cv` describes each one and its options.

* With `method = "twinning"`, only some proportions are allowed (for
  example, `2 / 3` but not `0.67`), and small strata are split at random. See
  `?balanced_initial_split` and `?balanced_initial_validation_split`.

* `balance_check()` reports standardized mean differences and the energy
  distance of balanced splits and of rsample's own splits and resamples. For
  objects that don't store `balance` or `exact` columns, such as those from
  `rsample::vfold_cv()`, name the columns to check. A three-way split is
  checked as the validation set and the testing set, each against the
  training set. With `variance_ratio = TRUE`, it also reports the variance
  ratio of each numeric variable. A variable with missing values also gets
  the SMD of its share missing, with `level` `"(missing)"`, and the energy
  distance uses all rows, with missing values encoded as with
  `na = "include"`. `summary()` gives one row per split, with
  its largest absolute SMD and its energy distance. `autoplot()` plots the
  result with ggplot2.

* Balanced folds work with `tune::fit_resamples()` and `tune::tune_grid()`,
  and balanced splits with `tune::last_fit()`. To tune on the validation set
  of a balanced three-way split, pass it to `rsample::validation_set()`
  first. To tune on balanced folds, attach balancesample with `library()`:
  rsample finds the function that made the folds by its name
  ([tidymodels/rsample#619](https://github.com/tidymodels/rsample/issues/619)).

* With tune 2.0.0 or later, a postprocessor that must be trained (for
  example, a tailor calibration) needs an inner split of the analysis set.
  * For balanced folds and a balanced training/testing split, the inner split
    uses the outer split's method and columns. If it can't be made (for
    example, if there are too few rows or groups), you get a warning and an
    empty calibration set, as in rsample.
  * Through `rsample::validation_set()`, rsample makes the inner split at
    random, so it is not balanced and can split groups
    ([tidymodels/rsample#620](https://github.com/tidymodels/rsample/issues/620)).
  * `tune::last_fit()` can't train such a postprocessor on any three-way
    split
    ([tidymodels/tune#1213](https://github.com/tidymodels/tune/issues/1213)).
