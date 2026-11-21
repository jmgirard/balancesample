# calibration falls back to an empty set when it can't split

    Code
      cal <- rsample::internal_calibration_split(folds$splits[[1]], args)
    Condition
      Warning:
      Cannot create calibration split; creating an empty calibration set.
      Caused by error in `balanced_initial_split()`:
      i In argument: `tidyselect::all_of(balance)`.
      Caused by error in `tidyselect::all_of()`:
      ! Can't subset elements that don't exist.
      x Element `not_a_column` doesn't exist.

# twinning calibration falls back to an empty set for a bad share

    Code
      cal <- rsample::internal_calibration_split(split, args)
    Condition
      Warning:
      Cannot create calibration split; creating an empty calibration set.
      Caused by error in `balanced_initial_split()`:
      ! `method = "twinning"` can only split off 1/r of the data, where r is a whole number.
      i Use `prop` = 1/2, 2/3, 3/4, 9/10, 1/4, or another value of 1 - 1/r or 1/r.

