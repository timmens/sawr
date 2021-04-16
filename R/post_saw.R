#' Post-SAW Estimation Procedure.
#'
#' Exports main function of the package, which performs the SAW procedure and
#' afterwards the Post-SAW estimation.
#'
#' @param y Matrix of labels. Has dimension T x N.
#' @param X List of feature matrices. The pth entry corresponds to the design
#' matrix of the pth covariate and has dimension T x N.
#' @param Z Instruments corresponding to the argument X. If NULL all X variables
#' are their own instrument.
#' @param time_effect Boolean indicating if a time effect is to be estimated.
#' @param id_effect Boolean indicating if an individual effect is to be
#' estimated and returned.
#' @param s_thresh Tuning parameter for the threshold lambda.
#' @export
fit_saw <- function(
  y, X, Z = NULL, time_effect = FALSE, id_effect = FALSE, s_thresh = NULL, return_info = FALSE
  ) {

  ## saw procedure
  saw_model <- saw_procedure(y, X, Z, s_thresh, return_info)

  x <- saw_model[["x"]]
  jump_locations <- saw_model[["jump_locations"]]
  gamma_hat <- saw_model[["gamma_hat"]]

  T <- nrow(y)
  N <- ncol(y)
  P <- ncol(x) / N
  M <- nrow(x)

  # update jump locations
  jump_locations <- set_entry_to_na_if_empty(jump_locations, P)
  jump_locations <- subtract_one(jump_locations)
  jump_locations <- drop_first_and_last_period(jump_locations)

  ## post-saw procedure
  post_data <- construct_data(y, x, jump_locations, time_effect)
  with_instrument <- !is.null(Z)
  if (with_instrument) {
    z <- saw_model[["z"]]
    instrument <- construct_data(y, z, jump_locations, time_effect)
    coeff <- internal_iv_reg(post_data$Y, post_data$X, instrument$X)
  } else {
    model <- stats::lm.fit(post_data$X, post_data$Y)
    coeff <- model$coefficients
  }

  coeff_list <- construct_coeff_list(jump_locations, coeff)
  beta_matrix <- construct_beta(coeff_list, jump_locations, M, coeff)

  ## return objects
  out <- list(
    beta_matrix = beta_matrix,
    jump_locations = jump_locations,
    coeff_list = coeff_list,
    gamma_hat = gamma_hat
  )

  if (time_effect) {
    out$time_effect <- estimate_time_effect(y, x, beta_matrix)
  }
  if (id_effect) {
    out$id_effect <- estimate_individual_effect(y, x, beta_matrix)
  }
  if (return_info) {
    out$info <- saw_model[["additional_information"]]
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
