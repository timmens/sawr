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

```R
optimal_model <- sawr::fit_saw_cv(y, X, n_folds=4, parallel=TRUE)
```


##### Prediction and instrumental variables

```R
X_new, X, Z, y <- get_endogeneous_data()

model <- sawr::fit_saw(y, X, Z)

pred <- predict(model, X_new)
```
