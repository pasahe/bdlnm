# bdlnm 0.1.0.10001

* Added support in `bdlnm()` for cross-basis interactions with fixed effects and basis-specific random effects.
* Added `combine_effects` to `bcrosspred()` to obtain predictions for fixed interaction levels and random-effect units. Interaction information is
  retained in the `basis_interaction` attribute.
* Added `select` to `plot.bcrosspred()` and `plot.optimal_exposure()` to select an interaction-specific result for plotting.
* Changed `optimal_exposure()` to operate primarily on predictions returned by `bcrosspred()`. Passing a fitted `bdlnm` model directly is deprecated.
* Changed `attributable()` to operate primarily on predictions returned by `bcrosspred()` and added support for fixed and random interaction-specific
  predictions. Passing a fitted `bdlnm` model directly is deprecated.
* Added the `columns` argument to `attributable()` for identifying exposure, outcome, date, filtering and interaction columns. The previous `name_*`
  arguments are deprecated, and the unused `range` argument was removed.
* Improved cross-basis coefficient-name handling and validation for models containing interactions and `f()` terms.
* Fixed `attributable()` so that, when a filter is supplied, total attributable numbers are calculated using all non-missing cases in the filtered period, including cases at time points where a daily attributable number is unavailable. Also corrected the validation for filters that leave no valid attributable numbers.
* Updated tests, documentation and vignettes for the prediction-first workflow and interaction models.

# bdlnm 0.1.0

- Initial CRAN release
