#' SAW Estimation Procedure
#'
#' @param formula Formula object. Intercept is eliminated by taking first
#' differences.
#' @param Z Instruments corresponding to X matrix in formula object. If NULL all
#' X variables are their own instrument.
#' @param time_effect Boolean indicating if the method needs to control for a
#' common time trend. This trend will then be estimated and returned.
#' @param id_effect Boolean indicating if an individual effect is to be
#' estimated and returned.
#' @param s.thresh Tuning parameter for the threshold lambda.
#' @export
saw_fun <- function(
  formula,
  Z = NULL,
  time_effect = FALSE,
  id_effect = FALSE,
  s.thresh = NULL
  ) {

  # T <- 33
  # N <- 300
  # beta <- make_beta(T, S=2)
  # data <- dgp3(T, N, beta$beta)
  # formula <- data$Y ~ data$X
  # s.thresh <- NULL

  ## saw procedure
  saw_model <- sawr:::BKSGL.pdm.default(formula, s.thresh)

  x.all.matrix <- saw_model$x.all.matrix
  y.matrix <- saw_model$y.matrix
  jump_locations <- saw_model$tausList

  T <- base::nrow(y.matrix)
  N <- base::ncol(y.matrix)
  P <- base::ncol(x.all.matrix) / N
  M <- base::nrow(x.all.matrix)

  # update jump locations
  jump_locations <- sawr:::set_entry_to_na_if_empty(jump_locations, P)
  jump_locations <- sawr:::subtract_one(jump_locations)
  jump_locations <- sawr:::drop_first_and_last_period(jump_locations)

  ## post-saw procedure
  data <- sawr:::construct_data(y.matrix, x.all.matrix, jump_locations, time_effect)
  if (is.null(Z)) {
    model <- stats::lm.fit(data$X, data$Y)
  } else {
    instrument <- sawr:::construct_data(y.matrix, Z, jump_locations, time_effect)
    model <- AER::ivreg(data$Y ~ data$X - 1, ~ instrument$X)
  }
  coeff <- model$coefficients

  coeff_list <- sawr:::construct_coeff_list(jump_locations, coeff)
  beta_matrix <- sawr:::construct_beta(coeff_list, jump_locations, M, coeff)

  out <- list(
    betaMat = beta_matrix,
    tausList = jump_locations,
    coeffList = coeff_list,
    X = data$X,
  )

  if (time_effect) {
    out["time_effect"] <- sawr:::estimate_time_effect(y.matrix, x.all.matrix, beta_matrix)
  }
  if (id_effect) {
    out["id_effect"] <- sawr:::estimate_individual_effect(y.matrix, x.all.matrix, beta_matrix)
  }

  return(out)
}

#' sawr
#'
#' This package contains the implementations of the panel data estimation
#' procedure suggested in Bada et al. (2020), titled 'A Wavelet Method for Panel
#' Models with Jump Discontinuities in the Parameters'
#' @references Bada et al. (2020). A Wavelet Method for Panel Models with Jump
#' Discontinuities in the Parameters'
#' @docType package
#' @name sawr
NULL
