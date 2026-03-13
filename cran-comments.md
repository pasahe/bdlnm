## Resubmission

This is a resubmission. The following issues raised by the CRAN manual reviewer have been addressed:

* Expanded the DESCRIPTION field to better explain the package functionality and implemented methods.

* Added the relevant references in the required format.

* Added `\value` tags to `plot.bcrosspred.Rd` and `plot.optimal_exposure.Rd` documenting no return value of those functions.

* Eliminated `:::` across all documentation examples by exporting the internal utility function `check_inla()` with `@keywords internal`.

* Ensured `par()` settings in `R/plot.bcrosspred.R` are restored using an immediate `on.exit()` call after being changed.

---

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

* Possibly misspelled words in DESCRIPTION: DLNM is a well-known acronym for Distributed Lag Non-linear Models and is not a misspelling.
  
* Cross-references to INLA functions in the documentation (e.g. [inla()]) cannot be checked when INLA is not installed, which is expected.

* The non-CRAN Suggested package INLA has been extensively tested with bdlnm locally and in GitHub Actions for Linux, Windows, and macOS. The needed repository specification is included in the package DESCRIPTION:

```
Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/stable
``` 