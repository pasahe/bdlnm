#' London temperature and mortality data
#'
#' Time series data containing the temperature and mortality from 1990 to 2012. Mortality data is stratified for <75 years and 75+ years age groups.
#'
#' @format ## `london`
#' A tibble with 8.279 rows and 7 columns:
#' \describe{
#'   \item{date}{Date}
#'   \item{year}{Year}
#'   \item{dow}{Day of the week}
#'   \item{tmean}{Temperature mean}
#'   \item{mort_00_74}{Mortality in the age group <75 years}
#'   \item{mort_75plus}{Mortality in the age group +75 years}
#'   \item{mortality}{All mortality}
#'   ...
#' }
#' @source <https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata>
"london"
