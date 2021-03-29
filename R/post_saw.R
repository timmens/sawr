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
#' @param time_effect Boolean indicating if the method needs to control for a
#' common time trend. This trend will then be estimated and returned.
#' @param id_effect Boolean indicating if an individual effect is to be
#' estimated and returned.
#' @param s.thresh Tuning parameter for the threshold lambda.
#' @export
saw_fun <- function(
  y, X, Z = NULL, time_effect = FALSE, id_effect = FALSE, s.thresh = NULL
  ) {

#   source("../simulation-saw-paper/src/R/dgp.R")
#   T <- 33
#   N <- 200
#   time_effect = TRUE
#   beta <- make_beta(T, S=2)
#
#   # time-effect
#   # data <- dgp5(T, N, beta$beta)
#   # Z <- NULL
#   # endogeniety
#   data <- dgp2(T, N, beta$beta)  # corr is .4 maybe increase this a little bit
#   Z <- data$Z
#   y <- data$Y
#   X <- data$X
#   s.thresh <- NULL


  ## create formula object since internal SAW method can only use formulas
  names(X) <- paste0("x", 1:length(X))
  data <- X
  data$y <- y
  formula <- paste0("y ~ ", paste(names(X), collapse = " + "))
  formula <- as.formula(formula, env=list2env(data))

  # this is not working..... (2SLS in SAW)
  # else {
  #   # instrument case; use 2SLS approach
  #   XX <- sapply(X, function(x) as.vector(x))
  #   ZZ <- sapply(Z, function(z) as.vector(z))
  #   projection_model <- lm.fit(ZZ, XX)
  #
  #   residuals <- as.matrix(projection_model$residuals)
  #   residuals <- split(residuals, rep(1:ncol(residuals), each=nrow(residuals)))
  #   names(residuals) <- paste0("x", 1:length(X))
  #   residuals <- lapply(residuals, function(x) matrix(x, ncol=ncol(y)))
  #
  #   data <- residuals
  #   data$y <- y
  #   formula <- paste0("y ~ ", paste(names(residuals), collapse = " + "))
  #   formula <- as.formula(formula, env=list2env(data))
  # }

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
  post_data <- sawr:::construct_data(y.matrix, x.all.matrix, jump_locations, time_effect)
  if (is.null(Z)) {
    model <- stats::lm.fit(post_data$X, post_data$Y)
    coeff <- model$coefficients
  } else {
    instrument <- sawr:::construct_data(y.matrix, Z[[1]], jump_locations, time_effect)
    coeff <- internal_iv_reg(post_data$Y, post_data$X, instrument$X)
  }

  coeff_list <- sawr:::construct_coeff_list(jump_locations, coeff)
  beta_matrix <- sawr:::construct_beta(coeff_list, jump_locations, M, coeff)

  out <- list(
    betaMat = beta_matrix,
    tausList = jump_locations,
    coeffList = coeff_list,
    X = data$X
  )

  if (time_effect) {
    out$time_effect <- sawr:::estimate_time_effect(y.matrix, x.all.matrix, beta_matrix)
  }
  if (id_effect) {
    out$id_effect <- sawr:::estimate_individual_effect(y.matrix, x.all.matrix, beta_matrix)
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
