#' SAW Estimation Procedure
#'
#' @param formula A formula object; Note that "-1" is unnecessary since the method takes first differences.
#' @param dot A boolean indicating wether a dot operation is necessary given the underlying model.
#' @param s.thresh A tuning parameter.
#' @export
saw_fun <- function(formula, dot = FALSE, s.thresh = NULL) {
  results <- BKSGL.pdm.default(formula, s.thresh)

  x.all.matrix <- results$x.all.matrix
  y.matrix     <- results$y.matrix
  tausList     <- results$tausList

  T <- nrow(y.matrix)
  tausList <- lapply(tausList, function(tauVec) tauVec[!tauVec %in% c(1, T)])
  tausList <- lapply(tausList, function(tauVec) {
                        ifelse(length(tauVec) == 0, NA, tauVec)
                    })

  linear_model_data <- construct_data_for_linear_model(y.matrix, x.all.matrix,
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

  list(betaMat = betaMat, tausList = tausList, coeffList = coeffList,
       X = linear_model_data$X)
}
