# error if no models nor coefficients are supplied

    Code
      bcrosspred(cb, at = temp)
    Condition
      Error in `bcrosspred()`:
      ! At least 'model' or 'coef' must be provided

# error if both model and coefficients are supplied

    Code
      bcrosspred(cb, mod = mod, coef = coef, at = temp)
    Condition
      Error in `bcrosspred()`:
      ! argument 2 matches multiple formal arguments

# error if another kind of model is supplied

    Code
      bcrosspred(cb, mod = mod_2, at = temp)
    Condition
      Error in `bcrosspred()`:
      ! argument 2 matches multiple formal arguments

# error if inla model doesn't have `control.config=TRUE`

    Code
      bcrosspred(cb, mod = mod_2, at = temp)
    Condition
      Error in `bcrosspred()`:
      ! argument 2 matches multiple formal arguments

# error if 'coef' is provided without a ' model.link'

    Code
      bcrosspred(cb, coef = coef, at = temp)
    Condition
      Error in `bcrosspred()`:
      ! 'model.link' has to be provided if a model is not supplied.

