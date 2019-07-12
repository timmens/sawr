########### wrapper function for BKS
# we can give here a formel and specify the form of the unobserved heterogeneity effects



BKSGL.pdm.default <-  function(formula,
                               s.thresh = NULL,
                               off.set = NULL) {
  # formula = Y~X1-1
  
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
  PF.obj	    <- FUN.PformulaBKSG(formula = formula)
  nc 		      <- PF.obj$nc
  nr 		      <- PF.obj$nr
  P_Over	    <- PF.obj$P_Over
  #intercept   <- PF.obj$is_intercept  # not used anywhere else ....
  dat.matrix	<- PF.obj$data.in.list
  x.all.matrix <- PF.obj$x.all.matrix # NEW (10.5,, for the third part of the estimation step one needs the data)
  dat.dim 	  <- c(nr, nc, P_Over)
  
  
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
                      list(x.all.matrix = x.all.matrix))
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
    
    FUN.blk <-
      function(lk) {
        # lk runs from 1 to (nrEinsBasis-1) by 2 #lk =2
        #print(lk)
        if (lk == 1) {
          H00 <- crossprod(x) # fehlt hier nicht /(nT)????????????????
          svdH00 <- eigen(H00)
          rm(H00)
          A00 <-
            svdH00[[2]] %*% (diag(svdH00[[1]] ^ {
              -0.5
            }, P)) %*% t(svdH00[[2]]) * sqrt(n * T)
          rm(svdH00)
          XI_lk <- x %*% A00
          #print(t(XI_lk)%*%XI_lk/(n*T))
          blk <- t(XI_lk) %*% y / (n * T)#;rm(XI_lk)
          Philk <- EinsBases[, 1, drop = FALSE] %x% A00
        }
        else{
          #lk = 2
          Einslk.n <- rep.int(EinsBases[, lk], n)
          Hlk <- crossprod(x[Einslk.n == 1, ])
          inv.Hlk <- solve(Hlk)
          rm(Hlk)
          
          Einslkp1.n <- rep.int(EinsBases[, (lk + 1)], n)
          Hlkp1 <- crossprod(x[Einslkp1.n == 1, ])
          inv.Hlkp1 <- solve(Hlkp1)
          rm(Hlkp1)
          
          svdinvHlkhalb <- eigen(inv.Hlk + inv.Hlkp1)
          invHlkhalb <-
            svdinvHlkhalb[[2]] %*% (diag(svdinvHlkhalb[[1]] ^ {
              -.5
            }, P)) %*% t(svdinvHlkhalb[[2]])
          rm(svdinvHlkhalb)
          
          Al2km1 <-  inv.Hlk %*% invHlkhalb * sqrt(n * T)
          rm(inv.Hlk)
          Al2k <- inv.Hlkp1 %*% invHlkhalb * sqrt(n * T)
          rm(inv.Hlkp1)
          
          #XI_lk <- rbind(x[Einslk.n ==1,]%*%Al2km1, -x[Einslkp1.n ==1,]%*%Al2k)
          XI_lk <- (x * Einslk.n) %*% Al2km1 - (x * Einslkp1.n) %*% Al2k
          #print(t(XI_lk)%*%XI_lk/(n*T))
          blk <-
            t(XI_lk[Einslk.n == 1 |
                      Einslkp1.n == 1 , ]) %*% y[Einslk.n == 1 |
                                                   Einslkp1.n == 1 , ] / (n * T)
          rm(XI_lk)
          Philk <-
            EinsBases[, lk, drop = FALSE] %x% Al2km1 - EinsBases[, (lk + 1), drop = FALSE] %x%
            Al2k
        }
        blkPhilk <- rbind(t(blk), Philk)
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
    b.coef <- Blk
    
    
    if (!is.null(s.thresh)) {
      thresh <- (((2 * log(T * P)) / (n * T ^ {
        1 / Kappa
      })) ^ {
        Kappa / 2
      }) * s.thresh
      # thresh <- (((2 * log(T)) / (n * T ^ {
      #   1 / Kappa
      # })) ^ {
      #   Kappa / 2
      # }) * s.thresh 
      b.coef[abs(b.coef) < thresh] = 0
    }
    else{
      naiv.resid <- y - rowSums(x * rep.int(gammaTilde, n))
      #naiv.var.resid <- min(var(naiv.resid), 
      #                      (quantile(naiv.resid, prob = .75) - quantile(naiv.resid, prob = .25)) / 1.34898)
      naiv.var.resid <- var(naiv.resid)
      thresh = sqrt(c(naiv.var.resid)) * ((2 * log(T * P)) / (n * T ^ {
        1 / Kappa
      })) ^ {
        Kappa / 2
      }
      #thresh = sqrt(c(naiv.var.resid)) * ((2 * log(T)) / (n * T ^ {
      #  1 / Kappa
      #})) ^ {
      #  Kappa / 2
      #}
      
      if (recalibrateThreshhold) {
        b.coef.plugIn    <-  b.coef
        b.coef.plugIn[abs(b.coef.plugIn) < thresh] = 0
        gammahat.plugIn  <- t(matrix(Philk %*% b.coef.plugIn, P, T))
        residuals.plugIn <- y - rowSums(x * rep.int(gammahat.plugIn, n))
        var.resid.plugIn <- var(residuals.plugIn)
        thresh           <- sqrt(c(var.resid.plugIn)) * ((2 * log(T * P))/(n*T^{1 / Kappa}))^{Kappa / 2}
      }
      b.coef[abs(b.coef) < thresh] = 0
    }
    
    gammahat  <- t(matrix(Philk %*% b.coef, P, T))
    residuals <- y - rowSums(x * rep.int(gammahat, n))
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
        tausList[[p]] = NA
      }
    }
    
    tausList
    
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


plot.bkgsl <-
  function(BKSGL_Obj, dates, plot = TRUE, ...) {
    
    SAW_gamma  = BKSGL_Obj$SAW_gamma
    gammaTilde = BKSGL_Obj$gammaTilde
    bTildeCoef = BKSGL_Obj$bTildeCoef
    W          = BKSGL_Obj$W
    thresh     = BKSGL_Obj$thresh
    residuals  = BKSGL_Obj$residuals
    EinsBases  = BKSGL_Obj$EinsBases
    TDiff      = BKSGL_Obj$DeltaDatDim[1]
    N          = BKSGL_Obj$DeltaDatDim[2]
    P_Over     = BKSGL_Obj$DeltaDatDim[3]
    P          = P_Over %/% 2
    tausList   = BKSGL_Obj$tausList
    c_Coef     = BKSGL_Obj$c_Coef
    
    SAW_betaShifted = rbind(rep(NA, P), SAW_gamma[, 1:P])
    SAW_betaUnShifted = rbind(SAW_gamma[, (P + 1):(2 * P)], rep(NA, P))
    
    SAW_beta = cbind(SAW_betaShifted, SAW_betaUnShifted)
    
    par(mfcol = c(1, 2))
    matplot(
      x = seq(1, (TDiff + 1)),
      SAW_beta,
      typ = "s",
      ylab = expression(gamma),
      lwd = 4,
      lty = 1,
      xlab = 't-index',
      xaxt = 'n'
    )
    
    #lines(B1, typ = "S", lwd =2, col = 'gray', lty = 1)
    for (p in 1:P) {
      #p =1
      abline(v = tausList[[p]], lty = 2)
      axis(1, tausList[[p]], tausList[[p]])
    }
    title('SAW_gamma_t')
    
    matplot(x = seq(1, (TDiff + 1)), c_Coef, xlab = 't-index')
    abline(h = c(thresh,-thresh))
    title('c_L- Coefficents')
  }


