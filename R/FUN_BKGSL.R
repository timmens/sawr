########### wrapper function for BKS
# we can give here a formel and specify the form of the unobserved heterogeneity effects



BKSGL.pdm.default <-  function(formula,
                               s.thresh = NULL,
                               off.set = NULL) {
  # formula = Y~X1+X2

  # check fomula

  if (!class(formula) == "formula") {
    stop(
      "\n Argument >>formula<< needs a formula-object like
         y~x1+... where the elements are matrices"
    )
  }

  # names

  names  <- names(model.frame(formula))

  # prepare for FUN.default (give the data in a large matrix: (y, x1, ...xP) where y, xp are  NT x1 Vectors  )
  PF.obj	     <- FUN.PformulaBKSG(formula = formula)
  nc 		       <- PF.obj$nc
  nr 		       <- PF.obj$nr
  P_Over	     <- PF.obj$P_Over
  #intercept   <- PF.obj$is_intercept  # not used anywhere else ....
  dat.matrix	 <- PF.obj$data.in.list
  x.all.matrix <- PF.obj$x.all.matrix # NEW (10.5, for the third part of the estimation step one needs the data)
  y.matrix     <- PF.obj$y.matrix # NEW (01.08, for the third part of the estimation step)
  dat.dim 	   <- c(nr, nc, P_Over)


  colnames(dat.matrix) <- c(names[1], paste(names[-1], '_Shifted'), paste(names[-1], '_UnShifted'))

  # Estimation results

  # OLD VERSION (10.5)
  #model.est <-
  #  FUN.BKSGL.pdm(dat.matrix = dat.matrix,
  #                dat.dim    = dat.dim,
  #                s.thresh   = s.thresh,
  #                ...)
  ##

  # NEW VERSION (Also returns x.all.matrix for the third estimation step)
  model.est <- append(FUN.BKSGL.pdm(
                          dat.matrix = dat.matrix,
                          dat.dim    = dat.dim,
                          s.thresh   = s.thresh),
                      list(x.all.matrix = x.all.matrix,
                           y.matrix = y.matrix))
  return(model.est)
}


FUN.BKSGL.pdm <-
  function(dat.matrix,
           dat.dim,
           s.thresh,
           recalibrateThreshhold = FALSE)
  {
    # data and dimension
    y 	<- dat.matrix[,  1, drop = FALSE]
    x 	<- dat.matrix[, -1, drop = FALSE]
    T 	<- dat.dim[1]
    n	  <- dat.dim[2]
    P	  <- dat.dim[3]# or ncol(x)

    WavDesgin <- function(T)
    {
      q <- log(T) / log(2)
      if (q != as.integer(q))
        stop(c("T should satisfies 2^q for an integer q."))
      EinsBases = matrix(1, T, 1)
      lastEinsBases = EinsBases[, ncol(EinsBases)]
      k = sum(lastEinsBases) / 2
      while (k >= 1)
      {
        K = T / k
        eins = EinsBases[1:k, 1]
        nbr.zero = k * K
        zeros = rep.int(0, nbr.zero)
        EinsBasesvec  <- c(rep(c(eins, zeros), (K - 1)), eins)
        newEinsBases  <- matrix(EinsBasesvec, T, K)
        EinsBases     <- cbind(EinsBases, newEinsBases)
        lastEinsBases <- EinsBases[, ncol(EinsBases)]
        k             <- sum(lastEinsBases) / 2
      }
      EinsBases
    }

    # construct father wavelates basis
    # print("construct EinsBases")
    EinsBases	  <- WavDesgin(T)
    nrEinsBasis	<- ncol(EinsBases)

    # plot(EinsBases[, 5])

    FUN.blk <- function(lk){ # lk runs from 1 to (nrEinsBasis-1) by 2 #lk =2
      #print(lk)
      if(lk == 1){
        H00 <- crossprod(x)/(n*T)# this is the Q11 (P_unders x P_unders) matrix in the paper
        svdH00 <- eigen(H00); rm(H00)
        A00 <- svdH00[[2]]%*%(diag(svdH00[[1]]^{-0.5}, P))%*%t(svdH00[[2]]); rm(svdH00) # This is the W11  (P_unders x P_unders) matrix in the paper
        XI_lk <- x%*%A00
        #print(t(XI_lk)%*%XI_lk/(n*T))
        blk <- t(XI_lk)%*%y/(n*T)#;rm(XI_lk)
        Philk <- EinsBases[,1, drop = FALSE]%x%A00
      }
      else{#lk = 2
        Einslk.n    <- rep.int(EinsBases[, lk], n)
        Hl2km1I     <- sqrt(T/sum(EinsBases[, lk]))
        Ql2km1      <- crossprod(x[Einslk.n ==1,])/(n*T)*Hl2km1I^2 # this is the Ql2k-1 (P_unders x P_unders) matrix in the paper
        inv.Ql2km1  <- solve(Ql2km1); rm(Ql2km1)
        
        Einslkp1.n <- rep.int(EinsBases[, (lk+1)], n)
        Hl2kI    <- sqrt(T/sum(EinsBases[, lk+1]))
        Ql2k1      <- crossprod(x[Einslkp1.n ==1,])/(n*T)*Hl2kI^2 # this is the Ql2k (P_unders x P_unders) matrix in the paper
        inv.Ql2k1  <- solve(Ql2k1); rm(Ql2k1)
        
        svdinvQlkhalb <- eigen((inv.Ql2km1 + inv.Ql2k1))
        invQlkhalb    <- svdinvQlkhalb[[2]]%*%(diag(svdinvQlkhalb[[1]]^{-.5}, P))%*%t(svdinvQlkhalb[[2]]); rm(svdinvQlkhalb)
        
        Al2km1 <-  inv.Ql2km1%*%invQlkhalb; rm(inv.Ql2km1)      # (P_unders x P_unders) matrix
        Al2k   <-  inv.Ql2k1%*%invQlkhalb; rm(inv.Ql2k1)  # (P_unders x P_unders) matrix
        
        #XI_lk <- rbind(x[Einslk.n ==1,]%*%Al2km1, -x[Einslkp1.n ==1,]%*%Al2k)
        XI_lk <- (x*Einslk.n)%*%Al2km1*Hl2km1I - (x*Einslkp1.n)%*%Al2k*Hl2kI
        #print(t(XI_lk)%*%XI_lk/(n*T))
        blk <- t(XI_lk[Einslk.n ==1 | Einslkp1.n ==1 ,])%*%y[Einslk.n ==1 | Einslkp1.n ==1 ,]/(n*T); rm(XI_lk)
        Philk <- EinsBases[,lk, drop = FALSE]%x%Al2km1*Hl2km1I - EinsBases[,(lk+1), drop = FALSE]%x%Al2k*Hl2kI # this is the W matrix in the paper
      }
      blkPhilk <- rbind(t(blk), Philk )
      #XI_lk
    }

    BlkPhilk <- FUN.blk(1)
    #t(FUN.blk(1))%*%FUN.blk(2)
    for (lk in seq.int(2, (nrEinsBasis - 1), by = 2)) {
      BlkPhilk <- cbind(BlkPhilk, FUN.blk(lk))
    }
    Philk <- BlkPhilk[-1,]
    Blk <- BlkPhilk[1,]
    gammaTilde <- t(matrix(Philk %*% Blk, P, T))

    Kappa = 1 - log(log(n * T)) / log(n * T)
    FUN_Thresh <- function(n, T, P, Kappa) {sqrt(P)*((2 * log(T * P)) / (n * T ^ {1 / Kappa})) ^ {Kappa / 2}}
    
    b.coef <- Blk


    if (!is.null(s.thresh)) {
      thresh <- FUN_Thresh(n, T, P, Kappa)* s.thresh
      b.coef[abs(b.coef) < thresh] = 0
    }
    else{
      rep_gammaTilde = gammaTilde
      for (i in 2:n) rep_gammaTilde = rbind(rep_gammaTilde, gammaTilde)
      naiv.resid <- y - rowSums(x*rep_gammaTilde)
      
      #naiv.resid <- y - rowSums(x * rep.int(gammaTilde, n))
      naiv.var.resid <- var(naiv.resid)
      thresh = sqrt(c(naiv.var.resid)) * FUN_Thresh(n, T, P, Kappa)


      if (recalibrateThreshhold) {
        b.coef.plugIn    <-  b.coef
        b.coef.plugIn[abs(b.coef.plugIn) < thresh] = 0
        gammahat.plugIn  <- t(matrix(Philk %*% b.coef.plugIn, P, T))
        residuals.plugIn <- y - rowSums(x * rep.int(gammahat.plugIn, n))
        var.resid.plugIn <- var(residuals.plugIn)
        thresh           <- sqrt(c(var.resid.plugIn)) * FUN_Thresh(n, T, P, Kappa)
      }
      b.coef[abs(b.coef) < thresh] = 0
    }

    gammahat  <- t(matrix(Philk %*% b.coef, P, T))
    
    rep_gammahat = gammahat
    for (i in 2:n) rep_gammahat = rbind(rep_gammahat, gammahat)
    
    residuals <- y - rowSums(x*rep_gammahat)
    var.resid <- var(residuals)

    # detect postSAW jumps

    TDiff      = dat.dim[1]
    N          = dat.dim[2]
    P_Over     = dat.dim[3]
    P          = P_Over %/% 2

    DiagScaledMatforL  = diag(sqrt(TDiff / 2), TDiff)
    Event              = seq.int(3, (TDiff + 1), by = 2)
    Oddt               = seq.int(2, TDiff, by = 2)

    PhiLk = DiagScaledMatforL[, (Event - 2)] - DiagScaledMatforL[, Oddt]

    c_Coef  = matrix(0, nrow = (TDiff + 1), ncol = P)
    tausList = vector(mode = 'list')
    for (p in 1:P) {
      #p = 1
      c_Coef[Event, p] = t(PhiLk) %*% gammaTilde[, p] / (TDiff)
      c_Coef[Oddt, p]  = t(PhiLk) %*% gammaTilde[, p + P] / (TDiff)
      tausList[[p]] = seq(1, (TDiff + 1))[abs(c_Coef[, p]) >= thresh]
      if (length(tausList[[p]]) > 0) {
        names(tausList[[p]]) = paste('tau_', p, '_', 1:length(tausList[[p]]), sep = '')
      }
      else{
        tausList[[p]] = NULL
      }
    }

    return(list(
      SAW_gamma   = gammahat,
      gammaTilde  = gammaTilde,
      bTildeCoef  = Blk,
      W           = Philk,
      thresh      = thresh,
      residuals   = residuals,
      EinsBases   = EinsBases,
      DeltaDatDim = dat.dim,
      c_Coef      = c_Coef,
      tausList    = tausList
    ))
  }


plot.bkgsl <- function(BKGSD_Obj, dates, plot = TRUE, ...){#BKGSD_Obj =result
  
  
  SAW_gamma  = BKGSD_Obj$SAW_gamma
  gammaTilde = BKGSD_Obj$gammaTilde
  bTildeCoef = BKGSD_Obj$bTildeCoef 
  W          = BKGSD_Obj$W 
  thresh     = BKGSD_Obj$thresh
  thresh2     = BKGSD_Obj$thresh2
  residuals  = BKGSD_Obj$residuals 
  EinsBases  = BKGSD_Obj$EinsBases
  TDiff      = BKGSD_Obj$DeltaDatDim[1]
  N          = BKGSD_Obj$DeltaDatDim[2]
  P_Over     = BKGSD_Obj$DeltaDatDim[3]
  P          = P_Over%/%2
  tausList   = BKGSD_Obj$tausList
  tausList_hut = BKGSD_Obj$tausList_hut
  c_Coef     = BKGSD_Obj$c_Coef
  c_Coef_hut = BKGSD_Obj$c_Coef_hut
  
  
  
  SAW_betaShifted = rbind(rep(NA, P), SAW_gamma[, 1:P, drop = FALSE])
  colnames(SAW_betaShifted) = paste('EvenBetaTilde_', 1:ncol(SAW_betaShifted), sep = '')
  SAW_betaUnShifted = rbind(SAW_gamma[, (P+1):(2*P), drop = FALSE], rep(NA, P))
  colnames(SAW_betaUnShifted) = paste('OddBetaTilde_', 1:ncol(SAW_betaUnShifted), sep = '')
  SAW_beta = cbind(SAW_betaShifted, SAW_betaUnShifted)
  
  
  gammaTildeShifted = rbind(rep(NA, P), gammaTilde[, 1:P, drop = FALSE])
  colnames(gammaTildeShifted) = paste('EvenBetaHut_', 1:ncol(gammaTildeShifted), sep = '')
  gammaTildeUnShifted = rbind(gammaTilde[, (P+1):(2*P), drop = FALSE], rep(NA, P))
  colnames(gammaTildeUnShifted) = paste('OddBetaHut_', 1:ncol(gammaTildeUnShifted), sep = '')
  bothgammatilde = cbind(gammaTildeShifted, gammaTildeUnShifted)
  
  par(mfcol = c(1, 2))
  matplot(x = seq(1, (TDiff+1)), SAW_beta, typ = "s", ylab = expression(gamma), lwd = 4, lty = 1, xlab = 't-index', xaxt = 'n', ylim = c(min(c(c(SAW_beta), c(bothgammatilde), -thresh), na.rm = TRUE) ,max(c(c(SAW_beta), c(bothgammatilde), thresh), na.rm = TRUE)))
  legend('bottomleft', legend = colnames(SAW_beta), lty = 1, col = seq(ncol(SAW_beta)), lwd = 4)
  
  
  matpoints(bothgammatilde, type = "o")
  
  for(p in 1:P){#p =1
    if (length(tausList[[p]]) > 0){
      abline(v = tausList[[p]], lty = 2)
      axis(1, tausList[[p]], tausList[[p]])}
  }
  title('SAW_gamma_t')
  
  Event              = seq.int(3, (TDiff+1), by=2)
  Oddt               = seq.int(2, TDiff, by=2)
  
  matplot(x = seq(1, (TDiff+1)), col= 'gray', c_Coef, xaxt = 'n', xlab = 't-index', ylim = c(min(c(c(c_Coef), -thresh), na.rm = TRUE) ,max(c(c(c_Coef), thresh), na.rm = TRUE)))
  matpoints(x = Event, c_Coef[Event,], col = "blue")
  matpoints(x = Oddt, c_Coef[Oddt,], col = "green")
  legend('bottomleft', legend = c('even', 'odd'), col = c('blue', 'green'), lty = 9, lwd = 4)
  abline(h = c(thresh, -thresh))
  title('c_L- Coefficents')  
  for(p in 1:P){#p =1
    if (length(tausList[[p]]) > 0){
      abline(v = tausList[[p]], lty = 2)
      axis(1, tausList[[p]], tausList[[p]])}
  }  
  
}