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

Data originally published by Vicedo-Cabrera (2019)
<doi:10.1097/EDE.0000000000000982>.

## References

Vicedo-Cabrera A.M., Sera F., Gasparrini A. (2019). Hands-on tutorial on
a modeling framework for projections of climate change impacts on
health. *Epidemiology*, 30(3), 321-329.
<doi:10.1097/EDE.0000000000000982>.
