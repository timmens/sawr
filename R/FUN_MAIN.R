#' SAW Estimation Procedure
#'
#' @param formula A formula object; Note that "-1" is unnecessary since the method takes first differences.
#' @param dot A boolean indicating wether a common time trend is to be eliminated.
#' @param s.thresh A tuning parameter further adjusting the threshold lambda.
#' @export
saw_fun <- function(formula, dot = FALSE, s.thresh = NULL) {

  # formula = data$Y ~ data$X1 + data$X2
  results <- sawr:::BKSGL.pdm.default(formula, s.thresh)

  x.all.matrix <- results$x.all.matrix
  y.matrix     <- results$y.matrix
  tausList     <- results$tausList

  T <- base::nrow(y.matrix)
  tausList <- base::lapply(tausList, function(tauVec) {
                        if (base::length(tauVec)) {
                          tauVec
                        } else {
                          NA
                        }
                    })

  tausList <- base::lapply(tausList, function(tau_vect) tau_vect - 1)
  tausList <- base::lapply(tausList, function(tauVec) tauVec[!tauVec %in% c(1, T)])

  # tauVec<-tausList[[2]]; tauVec[!tauVec %in% c(1, T)]

  linear_model_data <- sawr:::construct_data_for_linear_model(y.matrix, x.all.matrix,
                                                              tausList, dot)
  lm_fit_model      <- stats::lm.fit(linear_model_data$X, linear_model_data$Y)
  coeff             <- lm_fit_model$coefficients

  posit   <- base::cumsum(base::sapply(tausList, function(tau_vect) base::sum(!base::is.na(tau_vect)) + 1))
  posit   <- c(0, posit)
  coeffList <- list()
  for (i in base::seq_along(posit)[-1]) {
    coeffList[[i-1]] <- coeff[(posit[i-1] + 1):posit[i]]
  }

  betaMat  <- sawr:::construct_beta(coeffList, tausList, base::nrow(x.all.matrix))

  list(betaMat = betaMat, tausList = tausList, coeffList = coeffList,
       X = linear_model_data$X)
}



#' swar
#'
#' This package contains the implementations of the panel data estimation procedure suggested in Bada et al. (2020), titled 'A Wavelet Method for Panel Models with Jump Discontinuities in the Parameters'
#' @references Bada et al. (2020). A Wavelet Method for Panel Models with Jump Discontinuities in the Parameters'
#' @docType package
#' @name sawr
NULL
