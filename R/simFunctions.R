


randomEffects <- function( DES, slabs, mu, cs, intercept, tau, E,
                           glabs = NULL, wfill = NULL ){
  
  nf <- length(wfill)
  if( nf == 1 ){
    if( wfill == 0 )return(glabs)
  }
  
  ni <- length(slabs)
  st <- table(slabs)
  it <- as.numeric(names(st))
  gi <- rep(0, ni)
  wf <- c(1:ni)
  
  if( !is.null(wfill) ){
    wf   <- wfill
  }
  nf <- length(wf)
  
  if(DES == 'beta'){
    gi[wf] <- intercept + rnorm( nf, mu[wf], sqrt(tau) )
    return( gi )
  }
  
  if( DES == 'unstr' ){
    
    c0 <- sum( diag(cs)[it]*st )/ni/E # - crossprod(mu)/ni             # total var minus from beta
    gi[wf] <- intercept + mu[wf] + rnorm( nf, 0, sqrt( c0 + tau ) )
    return( gi )
  }
  
  if(DES == 'cov'){
    
    ik   <- as.matrix( expand.grid(slabs,slabs) )
    cvk  <- matrix( cs[ ik ], ni, ni )  + diag( tau, ni )  # individual by indiviaual
    mc   <- matrix(intercept, 1, ni)
    
    if( ni > 1 ){
      
      if( !is.null(wfill) & nf < ni ){
        
        wf   <- wfill
        wcon <- c(1:ni)[-wf]
        gi <- glabs
        #   gi[ wcon ] <- glabs[ wcon ]
        
        tmp <- .conditionalMVN( matrix(glabs, 1, ni), mc, cvk, cdex = wf )
        mc  <- tmp$mu
        cvk <- tmp$vr
        cvk <- (cvk + t(cvk))/2
      }
      
      gi[wf] <- rmvnormRcpp( 1, mc, cvk )
      
  #    gi[wf] <- as.vector(  .rMVN(1, mc, sigma = cvk)  ) # individual
      
    }else{
      
      gi <- suppressWarnings( rnorm(1, intercept, sqrt(cvk[1]) ) )
    }
  }
  gi
}

updateG <- function( gg, ge, slab ){
  
  # ge is gi rounded to the nearest integer
  
  if( length(ge) == 0 )return( list( gi = gg, snew = slab ) )
  
  gk   <- gg
  gk[ gk < 0 ] <- 0
  snew <- rep(slab, round(gk) )
  gnew <- gk
  
  # inherit pre-existing effects & species labels
  lastTab <- table(slab)
  thisTab <- table(snew)
  lastTab <- lastTab[ names(lastTab) %in% names(thisTab) ]
  new     <- thisTab[ names(lastTab) ] - lastTab
  # gnew    <- rep(0, sum(new))
  # snew    <- rep( as.numeric(names(new)), new )
  
  plus <- which(new >= 0)
  if( length(plus) > 0 ){
    snew    <- rep( as.numeric(names(new)[plus]), new[plus] )
    gnew    <- rep(0, sum(new[plus]))
  }
  
  minus <- which(new < 0 )
  if( length(minus) > 0 ){
    
    nminus <- -new[minus]
    
    for(i in 1:length(nminus)){
      wi <- which( slab == as.numeric(names(nminus)[i]) )
      ni <- min( c(nminus[i], length(wi)) )
      slab[ wi[1:ni] ] <-  0
    }
    wk <- which(slab > 0)
    slab <- slab[wk]
    gg <- gg[wk]
  }
  
  gupdate <- c(gg, gnew)
  supdate <- c(slab, snew)
  ord     <- order(supdate)
  gupdate <- gupdate[ord]
  supdate <- supdate[ord]
  
  list( gi = gupdate, snew = supdate )
}


randomSurvivors <- function( survSpec, supdate ){
  
  stab <- as.numeric( names(survSpec) )
  snew <- numeric(0)
  for(i in 1:length(survSpec)){
    
    ik <- sort( sample( which( supdate == stab[i] ), survSpec[i] ) )
    snew <- c(snew, ik)
  }
  names(snew) <- rep( stab, survSpec )
  snew
}

simLoop <- function(DES, ninit, K, S, nt, cc = S*10, z = 20, betaMat, cs, ss, tau,
                    q, xb, intercept = 2, survMat, tstart = nt/2, CONDDISP = T){
  
  # cs <- t(beta)%*%vx%*%beta     -or- C = S'VS
  # xb <- x%*%beta                -or- K = ES
  
  E <- nrow(betaMat)
  n <- matrix(ninit, K, S)             # site by species
  rownames(n) <- 1:K
  colnames(n) <- 1:S
  disp <- rep(0, S)                    # disperser pool
  nsum <- n*0                          # average over nt - tstart time steps
  
  if( DES == 'unstr' ){
    muSpec <- betaMat[ drop=F, sample(nrow(betaMat), 1), ]
  }
  
  
  slab <- rep( 1:S, each = ninit )     # species label
  n0   <- length(slab)
  slabel <- matrix( 0, K, 2*S*cc )   # sites by individuals (bigger than needed)
  rownames(slabel) <- rownames(x)
  slabel[,1:n0] <- matrix( slab, K, n0, byrow=T) # 
  glabel <- slabel*0
  
  fsurv <- as.numeric( colnames(survMat) )
  nsurv <- as.numeric( rownames(survMat) )
  survMat <- survMat*0
  
  nscores <- matrix( NA, nt, 3 )
  colnames(nscores) <- c('n', 'disp', 'ps' )

  ntot <- matrix(0, nt, S)
  
  pbar <- txtProgressBar(min=1, max=nt ,style=1)
  
  
  for(t in 1:nt){  # time iteration
    
    ns <- rowSums(n)         # no. individuals by site
    ntot[t,] <- colSums(n)   # time by species
    
    newCols  <- numeric(0)
    disperse <- numeric(0)
    
    for(k in 1:K){  # all sites
      
      sk <- n*0
      mk <- rep(1:S, n[k,])  
      nn <- length(mk)
      
      glast <- glabel[k, 1:nn]            # individual effect before
      
      slast <- slabel[k,1:nn]
      glast <- glast[ slast != 0 ]
      slast <- slast[ slast != 0 ]
      if( sum(slast) == 0 )next
      
      glabs <- glast
      
      if( t == 1 ){
        wfill <- glabs <- NULL
      }else{
        wfill <-  which( glast == 0 )
        if(length(wfill) == 0)wfill <- 0
      }
      
      mu <- xb[k,slast]
      if( DES == 'unstr' )mu <- muSpec[ slast ]
      
      gi <- randomEffects( DES, slabs = slast, mu = mu, cs, intercept, tau,
                           E = E, glabs = glabs, wfill = wfill )
      
      ge <- round( gi + rnorm( length(gi), 0, sqrt( ss ) ) )
      wf <- which( ge > 0 )
      
      gi <- gi[wf]
      if( sum(gi) > cc*2 ){
        gi <- gi/sum(gi)*cc*2
      }
      
      ssgg <- updateG( gg = gi, ge[wf], slab = slast[wf] )
      gupdate <- ssgg$gi
      supdate <- ssgg$snew
      
      wfill <- which(gupdate == 0)
      
      if( length(wfill) > 0 ){
        
        mu <- xb[k,supdate]
        if( DES == 'unstr' )mu <- muSpec[ supdate ]
        
        gupdate <- randomEffects( DES, slabs = supdate, mu = mu, 
                             cs, intercept, tau, E = E, glabs = gupdate, wfill = wfill )
        wf <- which(gupdate > 0)
        gupdate <- gupdate[wf]
        supdate <- supdate[wf]
      }
      
      ge <- gupdate + rnorm( length(gupdate), 0, sqrt( ss ) ) 
      
      if( length(supdate) > ncol(slabel) ){
        newCols <- matrix( 0, K, length(supdate) )
        newCols[ ,1:ncol(slabel)] <- slabel
        newCols[k,1:length(supdate)] <- supdate
        slabel <- newCols
        newCols[ ,1:ncol(glabel)] <- glabel
        newCols[k,1:length(gupdate)] <- gupdate
        glabel <- newCols
      }
      
      tg <- table( supdate )
      
      #  survival
      ps <-  1 - pnorm( (sum(tg) - cc)/z )  
      su <- cc + round( sum(tg)*ps )
      su <- min( c(su, sum(tg) ) )
      
      ig <- sample( length(supdate), su )     # survivors
      id <- sample( ig, round(q*length(ig)) ) # survive & disperse
      sk <- supdate[ ig ]
      
      ts <- table( sk )
      it <- as.numeric(names(ts))
      
      survSpec <- tg*0
      survSpec[ names(ts) ] <- ts/tg[names(ts)]
      
      notd <- ig[ !ig %in% id ]
      snew <- supdate[ notd ]
      gnew <- gupdate[ notd ]
      
      
      if( length(id) > 0 ){
        
        sdisp <- supdate[ id ]
        gdisp <- gupdate[ id ]
        
        if( min(sdisp) == 0 )stop('sdisp == 0')
        if(length(disperse) > 0){
          disperse <- rbind( disperse, cbind(sdisp, gdisp) )
        }else{
          disperse <- cbind( sdisp, gdisp )
        }
      }
      
      glabel[k, ] <- 0
      slabel[k, ] <- 0
      
      if( length(snew) > 0 ){
        slabel[k, 1:length(snew)] <- snew
        glabel[k, 1:length(gnew)] <- gnew
      }
      
      if(t > tstart){
        
        surv  <- findInterval( survSpec, nsurv, all.inside = T )
        ndens <- findInterval( tg, fsurv, all.inside = T )
        stab  <- table(surv, ndens)
        rownames(stab) <- nsurv[ as.numeric(rownames(stab)) ]
        colnames(stab) <- fsurv[ as.numeric(colnames(stab)) ]
        survMat[ rownames(stab), colnames(stab) ] <- survMat[ rownames(stab), colnames(stab) ] + stab
      }
    }   # end K loop
    
    if( length(disperse) > 0 ){
      
      loc <- sample(1:K, nrow(disperse), replace = T )
      il  <- unique(loc)
      for( ik in il ){
        drow <- disperse[ drop=F, loc == ik,  ]
        wk   <- which( slabel[ik, ] == 0 )[ 1:nrow(drow) ]
        slabel[ik, wk ] <- drow[,1]
        
        slabs <- slabel[ik, 1:max(wk) ]
        glabs <- glabel[ik, 1:max(wk) ]
        
        mu <- xb[k,slabs]
        if( DES == 'unstr' )mu <- muSpec[ slabs ]
        
        if(CONDDISP){  # cond distribution for destination
      
          glabel[ik, 1:max(wk) ] <- randomEffects( DES, slabs, mu = mu, 
                                                   cs, intercept, tau, E = E, glabs, wfill = wk )
        }else{         # unconditional
          
          glabel[ik, wk ] <- randomEffects( DES, slabs[wk], mu = mu[wk], 
                                            cs, intercept, tau, E = E) 
        }
      }
    }
    
    ws <- which( slabel > 0, arr.ind = T )
    
    if( length(ws) == 0 ){
      
      if( sum(survMat) == 0 ){
        survByDens <- NULL
      }else{
        sr <- range( which( rowSums(survMat) > 0 ) )
        sc <- range( which( colSums(survMat) > 0 ) )
        survByDens <- survMat[ sr[1]:sr[2], sc[1]:sc[2] ]
      }
      nmean <- nsum/(t - tstart)
      
      return( list( nmean = nmean, ntot = ntot[1:t,], nlast = n*0, survByDens = survByDens ) )
    }
    
    
    
    sk <- tapply( ws[,1]*0 + 1, list( site = ws[,1], species = slabel[ws] ), sum )
    sk[ is.na(sk) ] <- 0
    n <- n*0
    n[ rownames(sk), colnames(sk) ] <- sk
    
    if(t > tstart) nsum <- nsum + n
    
    dd <- 0
    if(length(disperse) > 0)dd <- round(nrow(disperse)/sum(n), 4 )
    
    
    nscores[t, ] <- c( round(sum(n)/K, 1), dd, round(ps, 3) )
    
    setTxtProgressBar(pbar,t)
    
  } # time loop
  
  
  
  
  sr <- range( which( rowSums(survMat) > 0 ) )
  sc <- range( which( colSums(survMat) > 0 ) )
  survByDens <- survMat[ sr[1]:sr[2], sc[1]:sc[2] ]
  
  nmean <- nsum/(nt - tstart)
  
  list( nmean = nmean, ntot = ntot, nlast = n, survByDens = survByDens,
        nscores = nscores )
}

diversity <- function( siteBySpec ){ # shannon diversity 
  
  pmat   <- sweep( siteBySpec, 1, rowSums(siteBySpec), '/' )
  pmat[ pmat == 0 ] <- NA
  bySite <- -rowSums( pmat*log(pmat), na.rm=T )
  
  byTot  <- colSums( round(siteBySpec) )     # at least one individual
  rich   <- length( which(byTot > 1) )
  byTot  <- byTot/sum(byTot)
  bt     <- byTot[ byTot > 0 ]
  byTot  <- -sum( bt*log(bt) )
  if( is.na(byTot) )byTot <- 0
  
  list( bySite = bySite, byTot = byTot, richness = rich )
}


plotTime <- function( timeBySpec ){
  
  time <- 1:nrow(timeBySpec)
  S    <- ncol(timeBySpec)
  
  plot(time, timeBySpec[,1], type='l', ylim=c(1, 2*ninit*K), 
       log = 'y', xlab='', ylab='' )
  for(s in 2:S){
    lines(time, timeBySpec[,s], col=s)
  }
}

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}


.shadeInterval <- function(xvalues,loHi,col='grey',PLOT  = TRUE, add  = TRUE,
                           xlab=' ',ylab=' ', xlim = NULL, ylim = NULL, 
                           LOG = FALSE, trans = .5){
  
  #draw shaded interval
  
  loHi <- as.matrix(loHi)
  
  tmp <- smooth.na(xvalues,loHi)
  
  xvalues <- tmp[,1]
  loHi <- tmp[,-1]
  
  xbound <- c(xvalues,rev(xvalues))
  ybound <- c(loHi[,1],rev(loHi[,2]))
  if(is.null(ylim))ylim <- range(as.numeric(loHi))
  if(is.null(xlim))xlim <- range(xvalues)
  
  if(!add){
    if(!LOG)plot(NULL, xlim = xlim, ylim=ylim, 
                 xlab=xlab, ylab=ylab)
    if(LOG)suppressWarnings( plot(NULL,  xlim = xlim, ylim=ylim, 
                                  xlab=xlab, ylab=ylab, log='y') )
  }
  
  
  if(PLOT)polygon(xbound,ybound, border=NA,col=.getColor(col, trans))
  
  invisible(cbind(xbound,ybound))
  
}

smooth.na <- function(x,y){   
  
  #remove missing values
  #x is the index
  #y is a matrix with rows indexed by x
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  
  wy <- which(!is.finite(y),arr.ind   = TRUE)
  if(length(wy) == 0)return(cbind(x,y))
  wy <- unique(wy[,1])
  ynew <- y[-wy,]
  xnew <- x[-wy]
  
  return(cbind(xnew,ynew))
}


.tnorm <- function(n,lo,hi,mu,sig, tiny=0){   
  
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  if(length(lo) != length(mu)){
    print(length(lo))
    print(length(mu))
    stop()
  }
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] + tiny
  z
}

.conditionalMVN <- function( xx, mu, sigma, cdex, p=ncol(mu) ){  
  # xx, mu are matrices
  # cdex conditional for these variables, must come last
  # gdex condition on these variables
  
  if(ncol(xx) != ncol(sigma))stop('ncol(xx) != ncol(sigma)')
  if(ncol(mu) != ncol(sigma))stop('ncol(mu) != ncol(sigma)')
  if(max(cdex) > ncol(mu))stop('max(cdex) > ncol(mu)')
  
  # new order
  ci   <- (1:p)[-cdex]
  new  <- c(ci, cdex)
  cnew <- match(cdex, new)
  pnew <- 1:(p - length(cnew))
  
  cond <- try( 
    condMVNRcpp(cnew-1, pnew-1, xx[,new, drop=F], 
                mu[,new, drop=F], sigma[new,new]), T)
  
  if( !inherits( cond,'try-error') ){
    return(cond)
    
  }else{
    
    sinv <- solve( sigma[drop=F, pnew, pnew] )
    p1   <- sigma[drop=F, cnew,pnew]%*%sinv
    mu1  <- mu[drop=F, ,cnew] + t( p1%*%t( xx[drop=F, ,pnew] - mu[drop=F, ,pnew] ) )
    vr1 <- solveRcpp( sinv[drop=F, cnew, cnew] )
    
    return( list( mu = mu1, vr = vr1 ) )
  }
}

.rMVN <- function (nn, mu, sigma = NULL, sinv = NULL){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.null(sigma)){
    m <- ncol(sigma)
  }else if(!is.null(sinv)){
    m <- ncol(sinv)
  }else{
    stop( '.rMNV requires either sigma or sinv' )
  }
  
  if(length(mu) > 1){
    if( !is.matrix(mu) ) mu <- matrix( mu, nn, length(mu) )  # mu is a vector of length m
    if( ncol(mu) == 1 & nn == 1 )  mu <- t(mu)
    if( length(mu) == m & nn > 1) mu <- matrix( mu, nn, length(mu), byrow=T )
  }
  
  if(is.null(sinv)){          # from sigma
    
    vv <- try(svd(sigma),T)
    
    if( inherits(vv,'try-error') ){
      ev <- eigen(sigma, symmetric = TRUE)
      rr <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    } else {
      rr <- vv$v %*% (t(vv$u) * sqrt(vv$d)) 
    }
    
  }else{     # from sinv
    
    L  <- chol(sinv)  
    rr <- backsolve(t(L), diag(m), upper.tri = F) 
  }
  ps <- matrix(rnorm(nn * m), nn) %*% rr
  ps + mu 
}


.tnormMVNmatrix <- function(avec, muvec, smat, 
                            lo=matrix(-10000,nrow(muvec),ncol(muvec)), 
                            hi=matrix(10000,nrow(muvec),ncol(muvec)),
                            whichSample = c(1:nrow(smat))){
  
  # lo, hi must be same dimensions as muvec,avec
  # each sample is a row
  
  lo[lo < -10000] <- -10000
  hi[hi > 10000]  <- 10000
  
  if(max(whichSample) > length(muvec))
    stop('\nwhichSample outside length(muvec)\n')
  
  whichSample <- sample(whichSample) # randomize order
  
  nd <- dim(avec)
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, 
                       lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) ) 
  r[,whichSample] <- a[,whichSample]
  r
}