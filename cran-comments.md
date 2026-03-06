## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

* Possibly misspelled words in DESCRIPTION: DLNM is a well-known acronym for
  Distributed Lag Non-linear Models and is not a misspelling.

* The non-CRAN Suggested package INLA has been extensively tested with bdlnm
  locally and in GitHub Actions for Linux, Windows, and macOS.
  The needed repository specification is included in the package DESCRIPTION:

```
Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/stable
``` 