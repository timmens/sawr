#' SAW Estimation Procedure
#' 
#' @param formula A formula object; Note that "-1" is unnecessary since the method takes first differences.
#' @param data A dataframe containing the variables. 
#' @param s.thresh A tuning parameter.
#' @export   
saw_fun <- function(formula, data, dot, s.thresh = NULL) {
  results <- BKSGL.pdm.default(formula, s.thresh)
  
  x.all.matrix <- results$x.all.matrix
  tausList     <- results$tausList
  
  linear_model_data <- construct_data_for_linear_model(data$Y, x.all.matrix, 
                                                       tausList, dot)
  lm_fit_model      <- lm.fit(linear_model_data$X, linear_model_data$Y)
  coeff             <- lm_fit_model$coefficients
  
  posit   <- cumsum(sapply(tausList, function(tau_vect) sum(!is.na(tau_vect)) + 1))
  posit   <- c(0, posit)
  coeffList <- list()
  for (i in seq_along(posit)[-1]) {
    coeffList[[i-1]] <- coeff[(posit[i-1] + 1):posit[i]]
  }
  
  tausList <- lapply(tausList, function(tau_vect) tau_vect - 1) 
  betaMat  <- construct_beta(coeffList, tausList, nrow(x.all.matrix))
  
  list(betaMat = betaMat, tausList = tausList, coeffList = coeffList)
}