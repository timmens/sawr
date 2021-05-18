# sawr 

This R-package contains the implementation of the SAW estimation procedure, introduced
in _A Wavelet Method for Panel Models with Jump Discontinuities in the Parameters_
(link to paper will appear soon).


## Installation

```R
library("devtools")
devtools::install_github("https://github.com/timmens/sawr")
```


## Usage

##### Standard case

```R
X, y <- get_data()

model <- sawr::fit_saw(y, X, time_effect = TRUE)
```

##### Cross validation

In case of distrust of the default threshold one can apply a standard two stage (coarse
and fine grid) cross-validation procedure.

```R
cv_model <- sawr::fit_saw_cv(y, X, n_folds=4, parallel=TRUE)
```

##### Threshold iteration

If cross-validation is too expensive on can usually improve upon the default threshold
by setting the threshold as the minimum of the absolute value of the true coefficients.
Since the true coefficients are not known we start with the default threshold and update
the threshold with the minimum absolute value of the coefficients until convergence is
reached.

```R
iterated_model <- sawr::fit_saw_iter(y, X, max_iter=100)
```

> Remark. Parallelization might only be available for unix users.


##### Prediction and instrumental variables

```R
X_new, X, Z, y <- get_endogeneous_data()

iv_model <- sawr::fit_saw(y, X, Z)

pred <- predict(iv_model, X_new)
```
