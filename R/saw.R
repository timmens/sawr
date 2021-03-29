## SAW procedure main script
#
# Includes main wrapper function which produces estimates \tilde{\gamma} in the
# paper.


BKSGL.pdm.default <-  function(formula, s.thresh = NULL, off.set = NULL) {
  # formula = Y~X1+X2

  # check formula

  if (!class(formula) == "formula") {
    stop(
      "\n Argument >>formula<< needs a formula-object like
         y~x1+... where the elements are matrices"
    )
  }

  # names

  names  <- names(model.frame(formula))

  # prepare for FUN.default (give the data in a large matrix: (y, x1, ...xP) where y, xp are  NT x1 Vectors  )
  PF.obj	     <- sawr:::FUN.PformulaBKSG(formula = formula)
  nc 		       <- PF.obj$nc
  nr 		       <- PF.obj$nr
  P_Over	     <- PF.obj$P_Over
  #intercept   <- PF.obj$is_intercept  # not used anywhere else ....
  dat.matrix	 <- PF.obj$data.in.list
  x.all.matrix <- PF.obj$x.all.matrix # NEW (10.5, for the third part of the estimation step one needs the data)
  y.matrix     <- PF.obj$y.matrix # NEW (01.08, for the third part of the estimation step)
  dat.dim 	   <- c(nr, nc, P_Over)


  colnames(dat.matrix) <- c(names[1], paste(names[-1], '_Shifted'), paste(names[-1], '_UnShifted'), 'delta_theta')

  # Estimation results

  # OLD VERSION (10.5)
  #model.est <-
  #  FUN.BKSGL.pdm(dat.matrix = dat.matrix,
  #                dat.dim    = dat.dim,
  #                s.thresh   = s.thresh,
  #                ...)
  ##

  # NEW VERSION (Also returns x.all.matrix for the third estimation step)
  model.est <- sawr:::FUN.BKSGL.pdm(
    dat.matrix = dat.matrix,
    dat.dim    = dat.dim,
    s.thresh   = s.thresh
    )
  out <- append(model.est, list(x.all.matrix = x.all.matrix, y.matrix = y.matrix))
  return(out)
}


FUN.BKSGL.pdm <- function(
  dat.matrix,
  dat.dim,
  s.thresh,
  recalibrateThreshhold = FALSE
  ) {
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
