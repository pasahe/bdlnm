#' Italian city temperature and mortality data
#'
#' Daily mean temperature and mortality in 87 Italian cities between 1 January 2011 and 31 December 2021.
#' Mortality is provided for five age groups and for all ages combined.
#' The dataset is based on the data used in Masselot and Gasparrini (2025).
#'
#' @md
#' @docType data
#' @keywords datasets
#' @usage data(italy)
#'
#' @format ## `italy`
#' A data frame with 349,566 rows and 13 columns, ordered by city and date.
#' Each row represents one city-day observation:
#' \describe{
#'   \item{date}{Observation date, of class `Date`.}
#'   \item{city_code}{City identifier.}
#'   \item{city_name}{City name.}
#'   \item{lon}{Longitude of the city location in decimal degrees (WGS 84).}
#'   \item{lat}{Latitude of the city location in decimal degrees (WGS 84).}
#'   \item{pop}{City population supplied with the source data.}
#'   \item{tmean}{Daily mean temperature.}
#'   \item{deaths}{Total daily deaths, calculated by summing the five age groups.}
#'   \item{deaths_0044}{Daily deaths among people aged 0--44 years.}
#'   \item{deaths_4564}{Daily deaths among people aged 45--64 years.}
#'   \item{deaths_6574}{Daily deaths among people aged 65--74 years.}
#'   \item{deaths_7584}{Daily deaths among people aged 75--84 years.}
#'   \item{deaths_85}{Daily deaths among people aged 85 years and older.}
#' }
#'
#' @source Data from the application in the paper: Masselot and Gasparrini (2025).
#'
#' @references
#' Masselot P. and Gasparrini A. (2025). A new modelling framework for
#' multi-location studies in environmental epidemiology.
#' _Statistical Methods in Medical Research_. <doi:10.1177/09622802241313284>.
#'
"italy"
