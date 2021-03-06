construct_data_for_linear_model <- function(y.matrix, x.all.matrix, tausList, dot) {
  
  tau_indicator <- function(tau_vector, j, T) {
    tau_vector <- c(1, tau_vector, T)
    
    f = ifelse(j != 1, function(t) ifelse(t > tau_vector[j] && t <= tau_vector[j + 1], 1, 0),
               function(t) ifelse(t >= tau_vector[j] && t <= tau_vector[j + 1], 1, 0))
    
    sapply(1:T, f)
  }
  
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
  
  dot <- ifelse(dot, function(X) X - apply(X, 1, mean), function(X) X)
  delta <- function(X) X[-1, ] - X[-nrow(X), ]
    
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
  regressor_list <- lapply(regressor_list, function(X) delta(dot(X)))
  regressor_list <- lapply(regressor_list, unname)
  names(regressor_list) <- regressor_names
  
  data_list <- append(list(Y = delta(dot(y.matrix))), regressor_list)
  data_list <- lapply(data_list, as.vector)
  Y         <- data_list[[1]]
  X         <- data_list[-1]
  X         <- do.call(cbind, X)
  
  return(list(Y = Y, X = X))
}

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