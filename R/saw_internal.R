#' @noRd
feature_list_to_matrix <- function(X) {
  #' Convert list of P (TxN) matrices to single matrix with shape T x (N * P).
  #' Data is stored as ["x_{t, 11}, x_{t, 12}, ..., x_{t, PN}"].
  x <- do.call(cbind, X)
  return(x)
}


#' @noRd
data_to_matrix <- function(delta_y, shifted_x, unshifted_x, P, N) {
  #' Construct data matrix where first column contains `delta_y` and other
  #' columns correspond to `\underscore{X}`. The reference can be found in the
  #' paper in section 2.1.

  data <- cbind(delta_y, shifted_x, unshifted_x)
  data <- sapply(
    1:(2*P+1), function(z, i) c(z[,seq((i-1)*N+1,i*N)]), z=data, simplify=TRUE
    )
  data <- cbind(data, 1)
  return(data)
}

#' @noRd
transform_input_data <- function(y, X, Z) {
  #' Transform input data which is either a T x N matrix (y) or a list of T x N
  #' matrices (X and Z) to data objects expected from the saw method. That is
  #' we construct `\delta y` from y and `\underscore{X}`, `\underscore{Z}` from
  #' X and Z, respectively. See section 2.1 in the paper.

  N <- ncol(y)
  T <- nrow(y)
  P <- length(X)

  x <- feature_list_to_matrix(X)

  delta_y <- y[-1, ] - y[-T, ]
  shifted_x <- x[-1, ]
  unshifted_x <- -x[-T, ]

  data <- data_to_matrix(delta_y, shifted_x, unshifted_x, P, N)

  if (!is.null(Z)) {
    z <- feature_list_to_matrix(Z)
    shifted_z <- z[-1, ]
    unshifted_z <- -z[-T, ]

    instrument <- data_to_matrix(delta_y, shifted_z, unshifted_z, P, N)
    instrument <- instrument[, -1]
  } else {
    instrument <- NULL
    z <- NULL
  }

  dimensions = list(N=N, T=T, P=P)

  out <- list(data=data, instrument=instrument, x=x, z=z, dimensions=dimensions)
  return(out)
}


#' @noRd
compute_threshold <- function(delta_y, x, z, N, TT, PP, s_thresh, gamma_tilde, kappa=NULL) {
  #' Compute threshold parameter.
  #'
  #' This function computes the threshold parameter used to shrink the `b`
  #' coefficients to zero. The parameter is found in the paper under the label
  #' `\lambda_{NT}`.

  preliminary_threshold <- compute_preliminary_threshold(N, TT, PP, kappa)

  if (is.numeric(s_thresh)) {

    multiplier <- s_thresh

  } else if (s_thresh == "residual") {

    rep_gamma_tilde <- repmat(gamma_tilde, N)

    with_instrument <- !is.null(z)
    if (with_instrument) {
      naive_residual <- delta_y - rowSums(z*rep_gamma_tilde)
    } else {
      naive_residual <- delta_y - rowSums(x*rep_gamma_tilde)

    }
    multiplier <- sd(naive_residual)

  } else if (s_thresh == "smalln") {

    multiplier <- log(N) * sqrt(sqrt(TT / N))

  } else {

    stop("s_thresh has to be a number of 'smalln' or 'residual'")

  }

  thresh <- multiplier * preliminary_threshold

  return(thresh)
}

#' @noRd
compute_preliminary_threshold <- function(N, TT, PP, kappa) {

  if (is.null(kappa)) {
    # use asymptotic argument if no kappa is specified
    kappa <- 1 - log(log(N * TT)) / log(N * TT)
  }

  threshold <- sqrt(PP) * ((2 * log(TT * PP)) / (N * TT ^ (1 / kappa))) ^ (kappa / 2)
  return(threshold)
}


#' @noRd
repmat <- function(mat, rows=1, cols=1) {
  #' Repeat matrix in the both axis `rows` and `cols` times.
  #'
  out <- kronecker(matrix(1, rows, cols), mat)
  return(out)
}
