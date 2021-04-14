## Post-SAW procedure internal auxiliary methods
#
# Exports auxiliary functions used in the post-SAW procedure.


#' @noRd
data_to_formula <- function(y, X) {
  names(X) <- paste0("x", 1:length(X))
  data <- X
  data$y <- y
  formula <- paste0("y ~ ", paste(names(X), collapse = " + "))
  formula <- as.formula(formula, env=list2env(data))
}


#' @noRd
tau_indicator <- function(tau_vector, j, T) {
  tau_vector <- c(1, tau_vector, T)

  f = ifelse(j != 1, function(t) ifelse(t > tau_vector[j] && t <= tau_vector[j + 1], 1, 0),
             function(t) ifelse(t >= tau_vector[j] && t <= tau_vector[j + 1], 1, 0))

  out <- sapply(1:T, f)
  return(out)
}

#' @noRd
construct_regressors_p <- function(tau_vector, X) {
  Sp <- length(tau_vector)
  T  <- nrow(X)
  N  <- ncol(X)

  if (any(is.na(tau_vector))) return(list(X))
  else{
    return(lapply(1:(Sp+1),
            function(j) X * matrix(tau_indicator(tau_vector, j, T), nrow=T, ncol=N)))
  }
}


#' @noRd
construct_data <- function(y.matrix, x.all.matrix, tausList, time_effect) {

  # dot and delta operator
  .dot <- ifelse(time_effect, function(X) X - apply(X, 1, mean), function(X) X)
  .delta <- function(X) X[-1, ] - X[-nrow(X), ]

  P     <- length(tausList)
  cols  <- rep(paste0("X", 1:P), each = (ncol(x.all.matrix) / P))
  count <- sapply(tausList,
                  function(tau_vect) {
                    ifelse(any(is.na(tau_vect)), 0, length(tau_vect))
                    })

  regressor_names <- vector("character", P + sum(count))
  index <- 0
  for (i in 1:P) {
    for (j in 1:(count[i] + 1)) {
      index <- index + 1
      regressor_names[index] <- paste0("X", i, j)
    }
  }

  X_list         <- lapply(split.data.frame(t(x.all.matrix), cols), t)
  regressor_list <- mapply(construct_regressors_p, tausList, X_list, SIMPLIFY = FALSE)
  regressor_list <- unlist(regressor_list, recursive = FALSE)
  regressor_list <- lapply(regressor_list, function(X) .delta(.dot(X)))
  regressor_list <- lapply(regressor_list, unname)
  names(regressor_list) <- regressor_names

  data_list <- append(list(Y = .delta(.dot(y.matrix))), regressor_list)
  data_list <- lapply(data_list, as.vector)
  Y         <- data_list[[1]]
  X         <- data_list[-1]
  X         <- do.call(cbind, X)

  return(list(Y = Y, X = X))
}

#' @noRd
construct_beta <- function(coeffList, tausList, T, coeff) {
  stopifnot(length(coeffList) == length(tausList))

  if (length(coeffList) == 0) {
    # do some stuff for beta whatever
    result_list <- matrix(coeff, nrow=T, ncol=1)
    rownames(result_list) <- rep("X11", T)
  } else {

    construct_beta_for_single_regressor <- function(coeff_vect, tau_vect) {
      stopifnot(!any(is.na(tau_vect)) || length(coeff_vect) == 1,
                any(is.na(tau_vect) || (length(coeff_vect) == length(tau_vect) + 1)))
      if (any(is.na(tau_vect))) return(rep(coeff_vect, T))

      tau_vect <- c(0, tau_vect, T)
      tau_rep  <- diff(tau_vect)

      rep(coeff_vect, times = tau_rep)
    }

    result_list <- mapply(construct_beta_for_single_regressor, coeffList, tausList)

  }

  return(result_list)
}

#' @noRd
set_entry_to_na_if_empty <- function(jump_locations, P) {
  updated <- list()
  for (i in 1:P) {
    vec <- tryCatch(
      expr = {vec <- jump_locations[[i]]},
      error = function(e) NA
    )
    updated[[i]] = vec
  }
  return(updated)
}


#' @noRd
subtract_one <- function(jump_locations) {
  updated <- base::lapply(jump_locations, function(vec) vec - 1)
  return(updated)
}


#' @noRd
drop_first_and_last_period <- function(jump_locations) {
  updated <- base::lapply(jump_locations, function(vec) vec[!vec %in% c(1, T)])
  return(updated)
}


#' @noRd
construct_coeff_list <- function(jump_locations, coeff) {
  ## what am i doing here?
  posit <- base::cumsum(base::sapply(jump_locations, function(vec) base::sum(!base::is.na(vec)) + 1))
  posit <- c(0, posit)
  coeff_list <- list()
  for (i in base::seq_along(posit)[-1]) {
    coeff_list[[i-1]] <- coeff[(posit[i-1] + 1):posit[i]]
  }
  return(coeff_list)
}


#' @noRd
estimate_mean <- function(y.matrix, x.all.matrix, betaMat) {
  if (ncol(betaMat) > 1) return("This does not yet work for multivariate beta!")
  mu <- mean(y.matrix)
  mu <- mu - (rowMeans(x.all.matrix) %*% betaMat) / length(betaMat)
  return(mu)
}


#' @noRd
estimate_time_effect <- function(y.matrix, x.all.matrix, betaMat) {
  if (ncol(betaMat) > 1) return("This does not yet work for multivariate beta!")
  .mean <- estimate_mean(y.matrix, x.all.matrix, betaMat)
  estimate <- rowMeans(y.matrix) - c(.mean)
  estimate <- estimate - rowMeans(x.all.matrix) * betaMat
  return(estimate)
}


#' @noRd
estimate_individual_effect <- function(y.matrix, x.all.matrix, betaMat) {
  if (ncol(betaMat) > 1) return("This does not yet work for multivariate beta!")
  .mean <- estimate_mean(y.matrix, x.all.matrix, betaMat)
  estimate <- colMeans(y.matrix) - c(.mean)
  estimate <- estimate - (t(x.all.matrix) %*% betaMat) / length(betaMat)
  return(estimate)
}


#' @noRd
internal_iv_reg <- function(y, X, Z) {
  # standard
  left <- solve(t(Z) %*% X)
  right <- t(Z) %*% y

  coefficients <- left %*% right
  return(coefficients)
}


#' @noRd
return_description <- function() {
  # this function beautifully shows why people will switch to a different programming language than R
  description = "The function fit_saw exports a list with arguments\n 1. beta_matrix\n2. jump_locations\n 3. coeff_list\n 4. gamma_hat\n 5. DESCRIPTION (this string)\n\n 2. represents the jump locations indexed as `tau_{p, i}`, i.e. the i-th jump\nlocation of the p-th covariate.\n\n3. represents the estimated beta coefficients for each time interval where\nthe coefficient is constant.\n\n1. represents the repeated coefficient matrix using the coefficients from 3.\n\n4. represents the first stage estimation of `gamma` (see paper for details)."
  return(description)
}
