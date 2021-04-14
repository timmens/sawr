#' @noRd
construct_wavelet_design <- function(TT) {

  q <- log(TT) / log(2)
  if (q != as.integer(q)) {
    stop("TT has to be a power of two; i.e. TT = 2**q, q in N")
  }

  wavelet_bases = matrix(1, TT, 1)
  last_wavelet_basis = wavelet_bases[, ncol(wavelet_bases)]
  k = sum(last_wavelet_basis) / 2

  while (k >= 1) {
    K = TT / k
    eins = wavelet_bases[1:k, 1]
    zeros = rep.int(0, k * K)

    basis_vector <- c(rep(c(eins, zeros), (K - 1)), eins)
    new_basis <- matrix(basis_vector, TT, K)
    wavelet_bases <- cbind(wavelet_bases, new_basis)
    last_wavelet_basis <- wavelet_bases[, ncol(wavelet_bases)]

    k <- sum(last_wavelet_basis) / 2
  }

  return(wavelet_bases)
}


#' @noRd
compute_w_matrix_and_b_tilde <- function(
  index, delta_y, x, z, wavelet_bases, TT, N, PP
) {
  #' Computes the `W_{l,k}` and `\tilde{b}_{l,k}` given index `{l,k}`.
  #'

  selector_minus <- rep.int(wavelet_bases[, index], N)
  selector <- rep.int(wavelet_bases[, index + 1], N)

  H <- construct_h_matrices(index, wavelet_bases, TT)
  Q <- construct_q_matrices(index, x, z, N, TT, selector_minus, selector, H)
  A <- construct_a_matrices(index, Q, PP)
  W <- construct_w_matrix(index, wavelet_bases, A, H)

  caligraphic_data <- construct_caligraphic_data(index, x, z, N, selector_minus, selector, A, H)
  b_tilde <- compute_b_tilde(index, delta_y, N, TT, wavelet_bases, caligraphic_data, selector_minus, selector)

  out <- list("W"=W, "b_tilde"=b_tilde)
  return(out)
}


#' @noRd
construct_h_matrices <- function(index, wavelet_bases, TT) {
  #' Construct `H_{l, 2k}, H_{l, 2k -1 }` matrices.
  #'
  #' The definition of these matrices is found in section 2.2 of the paper.

  out <- NULL
  if (index != 1) {

    # `H_{l, 2k}`
    H <- sqrt(TT / sum(wavelet_bases[, index + 1]))
    # `H_{l, 2k - 1}`
    H_minus <- sqrt(TT / sum(wavelet_bases[, index]))

    out <- list("H"=H, "H_minus"=H_minus)

  }

  return(out)
}


#' @noRd
construct_q_matrices <- function(index, x, z, N, TT, selector_minus, selector, H) {
  #' Construct `\underscore{Q}` matrices.
  #'
  #' These matrices can be found in section 2.3 of the paper.
  #'
  #' Q_minus corresponds to `\underscore{Q}_{l, 2k-1}` while Q, except for the
  #' case index == 1, corresponds to `\underscore{Q}_{l, 2k}`.
  #'
  #' For parameter `H` check function `construct_h_matrices`.

  with_instrument <- !is.null(z)

  if (index == 1) {

    # construct `\underscore{Q}_{1, 1}`
    if (with_instrument) {
      Q <- (t(x) %*% z) / (N * TT)
    } else {
      Q <- (t(x) %*% x) / (N * TT)
    }
    out <- Q

  } else {

    select_minus <- selector_minus == 1
    select <- selector == 1

    # construct `\underscore{Q}_{l, 2k - 1}` and `\underscore{Q}_{l, 2k}`
    if (with_instrument) {

      Q_minus <- ((t(x[select_minus, ]) %*% z[select_minus, ]) * H[["H_minus"]]^2) / (N * TT)
      Q <- ((t(x[select, ]) %*% z[select, ]) * H[["H"]]^2) / (N * TT)

    } else {

      Q_minus <- ((t(x[select_minus, ]) %*% x[select_minus, ]) * H[["H_minus"]]^2) / (N * TT)
      Q <- ((t(x[select, ]) %*% x[select, ])  * H[["H"]]^2) / (N * TT)

    }

    out <- list("Q"=Q, "Q_minus"=Q_minus)
  }

  return(out)
}


#' @noRd
construct_a_matrices <- function(index, Q, PP) {
  #' Construct `A` matrices.
  #'
  #' These matrices can be found in section 2.3 of the paper.
  #'
  #' A_minus corresponds to `A_{l, 2k - 1}` while A corresponds to `A_{l, 2k}`
  #' with the exception of index == 1 where A corresponds to `A_{1, 1}`.
  #'
  #' For parameter `Q` check function `construct_q_matrices`.


  if (index == 1) {

    Q_inv <- solve(Q)
    A <- robust_matrix_square_root(Q_inv)
    out <- A

  } else {

    Q_inv <- solve(Q[["Q"]])
    Q_minus_inv <- solve(Q[["Q_minus"]])

    summand <- Q_inv + Q_minus_inv
    to_multiply <- robust_matrix_square_root(summand)

    A_minus <- Q_minus_inv %*% to_multiply
    A <- Q_inv %*% to_multiply

    out <- list("A"=A, "A_minus"=A_minus)
  }

  return(out)
}


#' @noRd
construct_caligraphic_data <- function(index, x, z, N, selector_minus, selector, A, H) {
  with_instrument <- !is.null(z)

  if (index == 1) {

    x_caligraphic <- x %*% A

    if (with_instrument) {
      z_caligraphic <- z %*% A
    } else {
      z_caligraphic <- NULL
    }

  } else {

    x_minus <- x * selector_minus
    x <- x * selector

    x_caligraphic <- x_minus %*% (A[["A_minus"]] * H[["H_minus"]]) - x %*% (A[["A"]] * H[["H"]])

    if (with_instrument) {
      z_minus <- z * selector_minus
      z <- z * selector
      z_caligraphic <- z_minus %*% (A[["A_minus"]] * H[["H_minus"]]) - z %*% (A[["A"]] * H[["H"]])
    } else {
      z_caligraphic <- NULL
    }
  }

  out <- list(x=x_caligraphic, z=z_caligraphic)
  return(out)
}


#' @noRd
construct_w_matrix <- function(index, wavelet_bases, A, H) {

  if (index == 1) {

    W <- wavelet_bases[, index, drop=FALSE] %x% A

  } else {

    left <- wavelet_bases[, index, drop=FALSE] %x% A[["A_minus"]]*H[["H_minus"]]
    right <- wavelet_bases[,index+1, drop=FALSE] %x% A[["A"]]*H[["H"]]
    W <- left - right

  }

  return(W)
}


#' @noRd
compute_b_tilde <- function(index, delta_y, N, TT, wavelet_bases, caligraphic_data, selector_minus, selector) {

  with_instrument <- !is.null(caligraphic_data[["z"]])
  if (with_instrument) {
    regressor <- caligraphic_data[["z"]]
  } else {
    regressor <- caligraphic_data[["x"]]
  }

  if (index == 1) {

    b_tilde <- t(regressor) %*% delta_y / (N * TT)

  } else {

    select <- selector_minus == 1 | selector == 1
    b_tilde <- t(regressor[select, ]) %*% delta_y[select, ] / (N * TT)

  }

  return(b_tilde)
}


#' @noRd
drop_imaginary_part_if_zero <- function(x) {

  cols = ncol(x)
  rows = nrow(x)

  z <- zapsmall(x)
  if (all(Im(z) == 0)) {
    out <- as.numeric(z)
    out <- matrix(out, rows, cols)
  } else {
    out <- x
  }
  return(out)
}


#' @noRd
robust_matrix_square_root <- function(mat) {
  #' Compute the square root of `mat`.
  #'
  #' This method also works in the case where `mat` is non-symmetric but is
  #' not accurate in this case. For more information seek: https://tinyurl.com/97ek4f7u
  #'

  decom <- eigen(mat)

  eigen_vectors <- decom[[2]]
  eigen_values <- decom[[1]]
  eigen_values <- (eigen_values + abs(eigen_values)) / 2  # remove negative entries

  lambda <- diag(eigen_values**(0.5), ncol(mat))

  out <- eigen_vectors %*% lambda %*% t(decom[[2]])

  out <- drop_imaginary_part_if_zero(out)
  return(out)
}
