# London temperature and mortality data

The dataset includes observed daily mean temperature and total number of
deaths in London between 2000 and 2011. Mortality data is stratified for
\<75 years and 75+ years age groups. The dataset is based on the data
used in Vicedo-Cabrera et al. (2019).

## Usage

``` r
data(london)
```

## Format

### `london`

A tibble with 8.279 rows and 7 columns:

- time:

  Date index

- date:

  Date

- year:

  Year

- dow:

  Day of the week

- tmean:

  Temperature mean

- mort_00_74:

  Mortality in the age group \<75 years

- mort_75plus:

  Mortality in the age group +75 years

- mort:

  All mortality

## Source

Data originally published at
<https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata>

## References

Vicedo-Cabrera AM, Sera F, Armstrong B, Gasparrini A. A hands-on
tutorial on a modelling framework for projections of climate change
impacts on health. Epidemiology. 2019;30(3):321-329. DOI:
10.1097/EDE.0000000000000982. PMID: 30829832.
