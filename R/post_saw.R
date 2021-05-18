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
#' @param s_thresh Tuning parameter for the threshold lambda. Default to
#'   "residual". Can be either numeric or "residual" or "smalln".
#' @param return_info Return additional info on model fit.
#' @export
fit_saw <- function(
  y, X, Z = NULL, time_effect = TRUE, id_effect = TRUE, s_thresh = "residual", return_info = FALSE
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

  class(out) <- c("saw_model", "list")
  return(out)
}


#' Predict with a fitted saw model
#'
#' @param model fitted saw model output from function fit_saw
#' @param X newdata of same shape as the one used to fit model
#' @export
predict.saw_model <- function(model, X) {

  beta_matrix <- model[["beta_matrix"]]

  te <- as.numeric(model[["time_effect"]])
  ie <- as.numeric(model[["id_effect"]])

  prediction <- 0
  for (p in ncol(beta_matrix)) {
    prediction <- prediction + X[[p]] * beta_matrix[, p]
  }

  if (length(te) != 0) {
    prediction <- prediction + te
  }
  if (length(ie) != 0) {
    prediction <- t(t(prediction) + ie)
  }

  return(prediction)
}


#' Iterative threshold finder
#'
#' @description Solve for threshold using iterative procedure which sets the
#' threshold as the smallest absolute coefficient in iteration until convergence.
#' @param y Matrix of labels. Has dimension T x N.
#' @param X List of feature matrices. The pth entry corresponds to the design
#' matrix of the pth covariate and has dimension T x N.
#' @param Z Instruments corresponding to the argument X. If NULL all X variables
#' are their own instrument.
#' @param time_effect Boolean indicating if a time effect is to be estimated.
#' @param id_effect Boolean indicating if an individual effect is to be
#' estimated and returned.
#' @param max_iter Maximum number of iterations. Default 20.
#' @param numerical_threshold Threshold when to stop iteration. Default 0.001.
#' @param choose_min Boolean, if true use minimum of absolute coefficients else max.
#' @param return_info Return additional info on model fit.
#' @export
fit_saw_iter <- function(
  y,
  X,
  Z = NULL,
  time_effect = TRUE,
  id_effect = TRUE,
  max_iter=20,
  choose_min=TRUE,
  numerical_threshold=0.001,
  return_info = FALSE
  ) {
  func <- ifelse(choose_min, min, max)
  thresh <- "residual"
  thresh_updated <- NULL
  for (k in 1:max_iter) {
    updated <- fit_saw(y, X, Z, time_effect, id_effect, thresh, return_info)
    coeffs <- unlist(updated[["coeff_list"]])
    thresh_updated <- func(abs(coeffs))

    if (k > 1 && abs(thresh_updated - thresh) < numerical_threshold) {
      break
    }

    thresh <- thresh_updated
  }

  result <- fit_saw(y, X, Z, time_effect, id_effect, thresh, return_info)
  return(result)
}


#' Return cross-validated model
#'
#' @description Find optimal threshold parameter using k-fold cross validation
#' and return model fitted using this threshold.
#' @param y Matrix of labels. Has dimension T x N.
#' @param X List of feature matrices. The pth entry corresponds to the design
#' matrix of the pth covariate and has dimension T x N.
#' @param Z Instruments corresponding to the argument X. If NULL all X variables
#' are their own instrument.
#' @param time_effect Boolean indicating if a time effect is to be estimated.
#' @param n_folds Number of folds to use in cross validation step. Default 4. It
#' makes sense that n_folds is multiple of n_cores.
#' @param grid_size Number of s_thresh candidates.
#' @param prefer_sparsity Boolean indicating whether the cross-validation.
#' procedure should place a naive penalty on non-sparse solutions.
#' @param parallel Boolean indicating if code is parallelized over folds.
#' @param n_cores How many cores to use for paralellization. Default to number
#' of available logical cores. It makes sense that n_folds is multiple of
#' n_cores.
#' @param return_info Return additional info on model fit.
#' @export
fit_saw_cv <- function(
  y,
  X,
  Z = NULL,
  time_effect = FALSE,
  n_folds=4,
  grid_size=20,
  prefer_sparsity=TRUE,
  sparsity_level=0.5,
  max_threshold=10,
  parallel=FALSE,
  n_cores=NULL,
  return_info = FALSE
  ) {

  id_effect <- FALSE
  N <- dim(X[[1]])[2]
  folds <- caret::createFolds(1:N, k=n_folds)

  # partial out arguments
  func <- function(y, X, Z, s_thresh) {
    return(fit_saw(y, X, Z, time_effect=time_effect, id_effect=id_effect, s_thresh=s_thresh, return_info=TRUE))
  }

  # first: coarse sweep oriented on asymptotic choice
  result <- func(y, X, Z, "residual")
  s_thresh = result[["info"]][["thresh"]]
  s_thresh_list <- polyspace(s_thresh / 2, 20 * s_thresh, grid_size, degree=2)

  errors <- internal_cv(y, X, Z, func, s_thresh_list, folds, parallel, n_cores, prefer_sparsity)

  # second: detailed sweep
  candidates <- s_thresh_list[order(errors)[1:4]]
  minimum <- max(0, min(candidates) - sd(candidates))
  maximum <- max(candidates) + sd(candidates)
  s_thresh_list <- seq(minimum, maximum, length.out = grid_size)

  errors <- internal_cv(y, X, Z, func, s_thresh_list, folds, parallel, n_cores, prefer_sparsity)

  index <- which.min(errors)
  if (prefer_sparsity) {
    q <- quantile(ecdf(errors), sparsity_level)
    index <- index + sum(errors[index:grid_size] < q)
  }

  model <- fit_saw(y, X, Z, time_effect, id_effect, s_thresh_list[index], return_info)
  model$s_thresh <- c(minimum, s_thresh_list[index], maximum)
  return(model)
}


#' @noRd
internal_cv <- function(y, X, Z, func, s_thresh_list, folds, parallel, n_cores, penalty) {

  task <- function(k) {
    # task for fold k
    y_train <- y[, -folds[[k]]]
    y_test <- y[, folds[[k]]]
    X_train <- lapply(X, function(x) x[, -folds[[k]]])
    X_test <- lapply(X, function(x) x[, folds[[k]]])
    if (is.null(Z)) {
      Z_train <- NULL
      Z_test <- NULL
    } else {
      Z_train <- lapply(Z, function(z) z[, -folds[[k]]])
      Z_test <- lapply(Z, function(z) z[, folds[[k]]])
    }

    inner_errors <- c()
    for (s in seq_along(s_thresh_list)) {
      model <- func(y_train, X_train, Z_train, s_thresh_list[s])
      pred <- predict(model, X_test)
      jumps <- model$jump_locations
      penalty <- ifelse(penalty, sqrt(length(unlist(jumps))), 0)
      inner_errors[s] <- mean((y_test - pred)**2) + penalty
    }

    return(inner_errors)
  }

  if (parallel) {

    n_cores <- ifelse(is.null(n_cores), parallel::detectCores(), n_cores)
    errors <- parallel::mcmapply(task, seq_along(folds), mc.cores = n_cores)

  } else {
    errors <- sapply(seq_along(folds), task)
  }

  errors <- rowMeans(errors)
  return(errors)
}

#' @noRd
polyspace <- function(start, stop, len, degree=2) {
  f <- function(x) x ** (1 / degree)
  f_inv <- function(x) x ** degree

  .start <- .Machine$double.eps
  .stop <- stop - start

  space <- f_inv(seq(f(.start),f(.stop), length.out = len))
  space <- space + start
  return(space)
}
