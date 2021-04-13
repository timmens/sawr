saw_procedure <- function(y, X, Z=NULL, s_thresh=NULL, return_info=FALSE) {

  ## input data preparation

  transformed_data <- sawr:::transform_input_data(y, X, Z)

  data <- transformed_data[["data"]]
  instrument <- transformed_data[["instrument"]]

  dimensions <- transformed_data[["dimensions"]]
  N <- dimensions[["N"]]
  T <- dimensions[["T"]]
  P <- dimensions[["P"]]
  dimensions <- c(T, N, P)


  ## estimation

  model <- sawr:::saw_procedure_internal(data, instrument, dimensions, s_thresh, return_info)


  ## return objects

  out <- append(model, transformed_data)
  return(out)
}


saw_procedure_internal <- function(data, instrument, dimensions, s_thresh, return_info) {

  delta_y <- data[,  1, drop = FALSE]
  x <- data[, -1, drop = FALSE]  # this corresponds to `\underscore{X}_{it}` in the paper
  z <- instrument

  T <- dimensions[1]
  N	<- dimensions[2]
  P	<- dimensions[3]

  TT <- T - 1
  PP <- 2 * P + 1  # this corresponds to `\underscore{P}` in the paper


  ## construct father wavelet basis

  wavelet_bases <- sawr:::construct_wavelet_design(TT)


  ## compute `W` matrix and `tilde{b}`

  iterator <- c(1, seq.int(2, (ncol(wavelet_bases) - 1), by = 2))

  W <- c()
  b_tilde <- c()
  for (index in iterator) {

    tmp <- sawr:::compute_w_matrix_and_b_tilde(index, delta_y, x, z, wavelet_bases, TT, N, PP)
    W <- cbind(W, tmp[["W"]])
    b_tilde <- rbind(b_tilde, tmp[["b_tilde"]])

  }


  ## compute `\tilde{gamma}` from `W` and `\tilde{b}`

  gamma_tilde <- t(matrix(W %*% b_tilde, PP, TT))


  ## shrink `b` coefficients to zero

  thresh <- sawr:::compute_threshold(delta_y, x, N, TT, PP, s_thresh, gamma_tilde)
  # if needed the threshold can be recalibrated here

  b_hat <- b_tilde
  b_hat[abs(b_hat) < thresh] = 0


  ## compute `\hat{\gamma}` coefficients using `\hat{b}`

  gamma_hat <- W %*% b_hat
  gamma_hat <- t(matrix(gamma_hat, PP, TT))

  rep_gamma_hat <- sawr:::repmat(gamma_hat, rows=N)

  ## detect (postSAW) jump locations

  # the following corresponds to section 3.2 in the paper


  # TODO:
  # 1. Find out what DiagScaledMatforL means and give better name
  # 2. Find out what PhiLk means and give better name

  DiagScaledMatforL <- diag(sqrt(TT / 2), TT)

  even_time_period <- seq.int(3, (TT + 1), by = 2)
  odd_time_period <- seq.int(2, TT, by = 2)

  PhiLk <- DiagScaledMatforL[, (even_time_period - 2)] - DiagScaledMatforL[, odd_time_period]

  c_tilde <- matrix(0, nrow = (TT + 1), ncol = P)

  jump_locations <- list()
  for (p in 1:P) {

    c_tilde[even_time_period, p] <- t(PhiLk) %*% gamma_tilde[, p] / (TT)
    c_tilde[odd_time_period, p] <- t(PhiLk) %*% gamma_tilde[, p + P] / (TT)

    c_hat <- abs(c_tilde[, p]) > thresh

    jump_locations[[p]] = seq(1, (TT + 1))[c_hat]  # how does this work?

  }

  # rename jump locations according to feature label
  for (p in 1:P) {

    if (length(jump_locations[[p]]) > 0) {
      names(jump_locations[[p]]) = paste('tau_', p, '_', 1:length(jump_locations[[p]]), sep = '')
    } else {
      jump_locations[[p]] = NULL
    }

  }


  ## compute additional quantities
  residuals <- delta_y - rowSums(x * rep_gamma_hat)


  ## return objects

  out <- list(
    jump_locations = jump_locations,
    gamma_hat = gamma_hat
  )

  if (return_info) {
    additional_information <- list(
      gamma_tilde = gamma_tilde,
      b_tilde = b_tilde,
      W = W,
      thresh = thresh,
      residuals = residuals,
      wavelet_bases = wavelet_bases,
      DeltaDatDim = dimensions,
      c_tilde = c_tilde
    )
    out <- append(out, additional_information)
  }

  return(out)
}
