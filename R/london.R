#' London temperature and mortality data
#'
#' The dataset includes observed daily mean temperature and total number of deaths in London between 2000 and 2012. Mortality data is stratified for <75 years and 75+ years age groups.
#'
#' @docType data
#' @keywords datasets
#' @usage data(london)
#'
#' @format ## `london`
#' A tibble with 8.279 rows and 7 columns:
#' \describe{
#'   \item{time}{Date index}
#'   \item{date}{Date}
#'   \item{year}{Year}
#'   \item{dow}{Day of the week}
#'   \item{tmean}{Temperature mean}
#'   \item{mort_00_74}{Mortality in the age group <75 years}
#'   \item{mort_75plus}{Mortality in the age group +75 years}
#'   \item{mort}{All mortality}
#' }
#'
#' @source <https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata
#'
#' @references
#' Vicedo-Cabrera AM, Sera F, Armstrong B, Gasparrini A. A hands-on tutorial on a modelling framework for projections of climate change impacts on health. Epidemiology. 2019;30(3):321-329. DOI: 10.1097/EDE.0000000000000982. PMID: 30829832.
#'
"london"
