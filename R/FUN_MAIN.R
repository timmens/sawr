#' SAW Estimation Procedure
#'
#' @param formula A formula object; Note that "-1" is unnecessary since the method takes first differences.
#' @param dot A boolean indicating wether a common time trend is to be eliminated.
#' @param s.thresh A tuning parameter.
#' @param ridge A stabilization parameter to avoid inverting (near) singular design matrices (t(X)X)
#' @export
saw_fun <- function(formula, dot=FALSE, s.thresh=NULL, ridge=10e-5) {
  results <- BKSGL.pdm.default(formula, s.thresh, ridge=ridge)

  x.all.matrix <- results$x.all.matrix
  y.matrix     <- results$y.matrix
  tausList     <- results$tausList

  T <- nrow(y.matrix)
  tausList <- lapply(tausList, function(tauVec) {
                        if (length(tauVec)) {
                          tauVec
                        } else {
                          NA
                        }
                    })

  tausList <- lapply(tausList, function(tau_vect) tau_vect - 1)
  tausList <- lapply(tausList, function(tauVec) tauVec[!tauVec %in% c(1, T)])

  linear_model_data <- construct_data_for_linear_model(y.matrix, x.all.matrix,
                                                       tausList, dot)
  X <- linear_model_data$X
  Y <- linear_model_data$Y

  coeff <- solve(t(X) %*% X) %*% t(X) %*% Y

  posit   <- cumsum(sapply(tausList, function(tau_vect) sum(!is.na(tau_vect)) + 1))
  posit   <- c(0, posit)
  coeffList <- list()
  for (i in seq_along(posit)[-1]) {
    coeffList[[i-1]] <- coeff[(posit[i-1] + 1):posit[i]]
  }

  betaMat  <- construct_beta(coeffList, tausList, nrow(x.all.matrix))

  list(betaMat = betaMat, tausList = tausList, coeffList = coeffList,
       X = linear_model_data$X)
}
