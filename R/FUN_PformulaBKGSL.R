## #=============================================================================
## FUN.Pformula calls FUN.with.trans
##
## Takes: fomula = a formula-object
##        effect = character-object ("none", "time", "individual", "twoways")
##
## Gives: A List with P+1 components
##        (first for depend. variable second for indep. variables):
##        Each Component is a List with 6 Components:
##        1. Tr  (one of "none", "time", "individual", "twoways")
##        2. I   (TRUE, if intercept, FALSE, if not)
##        3. ODM (Original Data-Matrix)
##        4. TDM (Transformed Data-Matrix)
##        5. TDV (TDM as a vector)
##        6. TRm is a list with 3 Components:
##           1. OVc  (Overall constant)
##           2. InC  (Individual constants)
##           3. TiVC (Time variaing constants)
## #==============================================================================

FUN.PformulaBKSG <- function(formula, timeEffect = FALSE) {
  data.fra <- model.frame(formula)
  dat.term <- attr(data.fra, "terms")

  ## Construct data from formula
  ## dim(response) == TxN == dim(y.matrix) == TxN:
  y.matrix <- model.response(data.fra, "numeric")

  ## Check Data: Each Variable has to be a matrix
  if(!is.matrix(y.matrix)) stop("Each Variable has to be a TxN-matrix.")
  N  <- ncol(y.matrix)
  T  <- nrow(y.matrix)

  ## 1)Extract Regressors
  ## 2)Check the presence of a 'intercept' in the formula
  ## 3)And built the x.all.matrix (Tx[N*P])

  regressors.mat <- model.matrix(dat.term, data.fra)
  is.intercept   <- (colnames(regressors.mat)[1] == "(Intercept)")
  if(is.intercept) {
    x.all.matrix <- regressors.mat[,-1]
  } else {
    x.all.matrix <- regressors.mat
  }

  ## Dimension parameters
  NT <- N * T
  P  <- as.integer(ncol(x.all.matrix) / N)
  if(P != ncol(x.all.matrix) / N) stop("All the regressors must have the same dimension as the response variable Y")
  ### ^ Is that implication necessarily correct?

  Delta_y.matrix         <-   y.matrix[-1, ] - y.matrix[-T, ]
  Shifted_x.all.matrix   <-   x.all.matrix[-1, ]
  UnShifted_x.all.matrix <-  -x.all.matrix[-T, ]

  data.all.mat  <- cbind(Delta_y.matrix, Shifted_x.all.matrix, UnShifted_x.all.matrix)

  ## Write the response variable, Y, and the augmented 'p' regressors, X, in a list,
  ## where each component contains one of 2*p+1 TN-vector
  data.in.list <- sapply(1:(2*P+1), function(z, i) c(z[,seq((i-1)*N+1,i*N)]), z = data.all.mat, simplify = TRUE)

  model.in.list <- list(
    data.in.list = data.in.list,
    nr           = T - 1,
    nc           = N,
    P_Over       = 2*P,
    is_intercept = is.intercept,
    x.all.matrix = x.all.matrix,
    y.matrix     = y.matrix
    )

  model.in.list
}
