# sawr 

This R-package contains the implementation of the SAW estimation procedure, introduced
in _A Wavelet Method for Panel Models with Jump Discontinuities in the Parameters_
(link to paper will appear soon).


## Installation

```{R}
library("devtools")
devtools::install_github("https://github.com/timmens/sawr")
```


## Usage

```{R}
X, y <- get_data()

model <- sawr::fit_saw(y, X, time_effect = TRUE)
```
