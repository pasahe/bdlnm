# error if no models nor coefficients are supplied

    Code
      bcrosspred(bdlnm:::cb_london, at = temp)
    Condition
      Error in `bcrosspred()`:
      ! At least 'model' or 'coef' must be provided

# error if both model and coefficients are supplied

    Code
      bcrosspred(bdlnm:::cb_london, mod = bdlnm:::mod_london, coef = bdlnm:::coef_london,
      at = temp)
    Condition
      Error in `bcrosspred()`:
      ! argument 2 matches multiple formal arguments

# error if another kind of model is supplied

    Code
      bcrosspred(bdlnm:::cb_london, mod = mod_2, at = temp)
    Condition
      Error in `bcrosspred()`:
      ! argument 2 matches multiple formal arguments

# error if 'coef' is provided without a ' model.link'

    Code
      bcrosspred(bdlnm:::cb_london, coef = bdlnm:::coef_london, at = temp)
    Condition
      Error in `bcrosspred()`:
      ! 'model.link' has to be provided if a model is not supplied.

