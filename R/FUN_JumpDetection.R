DetectPostSAWJumps = function(BKGSD_Obj, plotCoef = TRUE){
  SAW_gamma  = BKGSD_Obj$SAW_gamma
  gammaTilde = BKGSD_Obj$gammaTilde
  bTildeCoef = BKGSD_Obj$bTildeCoef 
  W          = BKGSD_Obj$W 
  thresh     = BKGSD_Obj$thresh
  residuals  = BKGSD_Obj$residuals 
  EinsBases  = BKGSD_Obj$EinsBases
  TDiff      = BKGSD_Obj$DeltaDatDim[1]
  N          = BKGSD_Obj$DeltaDatDim[2]
  P_Over     = BKGSD_Obj$DeltaDatDim[3]
  P          = P_Over%/%2
  
  
  
  DiagScaledMatforL = diag(sqrt(TDiff/2), TDiff)
  Event              = seq.int(3, (TDiff+1), by=2)
  Oddt               = seq.int(2, TDiff, by=2)
  
  PhiLk = DiagScaledMatforL[, (Event-2)] - DiagScaledMatforL[, Oddt]
  
  c_Coef  = matrix(0, nrow = (TDiff+1), ncol = P)
  tausList = vector(mode = 'list')
  for(p in 1:P){#p = 1
    c_Coef[Event, p] = t(PhiLk)%*%gammaTilde[, p]/(TDiff)
    c_Coef[Oddt, p]  = t(PhiLk)%*%gammaTilde[, p+P]/(TDiff)
    tausList[[p]] = seq(1,(TDiff+1))[abs(c_Coef[, p])>= thresh]
    if(length(tausList[[p]]) > 0){
      names(tausList[[p]]) = paste('tau_',p, '_', 1:length(tausList[[p]]), sep = '')
    }
  }

  if(plotCoef){
    matplot(c_Coef)
    abline(h = c(thresh, -thresh))
  }
  tausList
}







 