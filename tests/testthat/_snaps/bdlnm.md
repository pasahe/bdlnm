# bdlnm errors when required arguments are missing or inappropiate

    A `data` data.frame must be provided.

---

    A basis of class <'crossbasis'> or <'onebasis'> must be provided to `basis`.

---

    Failed to fit INLA model: *** Fail to get good enough initial values. Maybe it is due to something else.
    i Check the model formula, family, data, and additional options passed via `...`.

---

    `sample.arg` must be a "list".

---

    Failed to draw posterior samples via `INLA::inla.posterior.sample()`: unused argument (not_an_argument = 5)
    i Check the additional options passed via `sample.arg`.

---

    Failed to fit INLA model: unused argument (not_an_argument = 5)
    i Check the model formula, family, data, and additional options passed via `...`.

# bdlnm throws informative error if control.compute config = FALSE is passed

    `control.compute(config = FALSE)` cannot be provided: `config` must be "TRUE" to use `inla.posterior.sample()`.

