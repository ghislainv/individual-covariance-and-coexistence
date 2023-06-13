## Load libraries
library("here")
library('Rcpp')
library('RcppArmadillo')

## Source code
Rcpp::sourceCpp(here("src", "cppFns.cpp"))
source(here("R", "simFunctions.R"))

################## a few notes on dimensionality
E <- 100
S <- 5
n <- 100
x <- matrix( rnorm(n*E), n, E )
x <- sweep(x, 2, colMeans(x), '-' )
#x <- sweep(x, 2, apply(x, 2, sd), '/' )
v <- cov(x)
s <- matrix( rnorm(E*S), E, S)
s <- sweep(s, 1, rowMeans(s), '-' )
cs <- cov(s)

t(s)%*%v%*%s

#v <- diag(E)
cs <- t(s)%*%v%*%s
crossprod(s)
E*cov(s)

##################

K     <- 500     
ninit <- 2
E <- 1         # dimensionality E (columns in X (E in ms))
S <- 100          # no. species
varBeta <- 1    # determines magnitude of variation in XB
varX    <- 1

intercept <- 1.1

cc <- 200
z  <- cc/2         #  standard deviation in survival probit
q  <- .1          # dispersal fraction
ss <- .001         # residual variance sigma^2 in ms
tau <- 1e-2        # nugget variance 


nt     <- 100
nsim   <- 10
nstart <- nt/2


beta <- matrix( rnorm(E*S, 0, sqrt(varBeta) ), E, S )     # X, B, V (matrix S in ms)
vx   <- cov( rmvnormRcpp(10, rep(0, E), varX*diag(E) ) )  
x    <- rmvnormRcpp(K, rep(0, E), vx )                    # matrix E in ms
x    <- sweep( x, 2, colMeans(x), '-' )                   # centered, not standardized
vx   <- cov(x)                                            # matrix V in ms
cs <- t(beta)%*%vx%*%beta                                 # C = S'VS
xb <- x%*%beta                                            # XB (ES in ms)
colnames(x) <- rownames(beta) <- paste('q',1:nrow(beta), sep='')

words <- paste('There are ', ninit, ' initial individuals for each of $S$ = ', S, ' species on $K$ = ', K, 
               ' sites. Parameter values are $E$ = ', E, ', $c$ = ', cc, ', $q$ = ' , q, 
               ', $z$ = ', z, ', $sigma^2$ = ', ss, ', $tau^2$ = ', tau,
               ', Var[vec($\bB$)] = ', varBeta, '$\bI_{ES}$, Var[vec($\bX$)] = ', varX, '$\bI_{nE}$', sep='')
words

#par(mfrow=c(2,3), bty='n', mar=c(4,4,1,1), oma=c(4,3,2,4))



DESIGN <- 'beta'    # XB (full knowledge)

DESIGN <- 'cov'     # spp differences summarized by C

DESIGN <- 'unstr'   # noisy individuals

n <- matrix(ninit, K, S)
ntot <- matrix(0, nt, S)

colnames(beta) <- colnames(n) <- colnames(ntot) <- paste('s',1:S, sep='')
rownames(x) <- rownames(n) <- paste('site',1:K,sep='_')


# three designs for 5 levels of E

eseq <- round( c(4, 8, 12, 24, 32, 44, 55, 65, 80, 100, 120, E) )

eseq <- round( c(4, 8, 12, 120,  E) )

eseq <- c(seq(1, 10, by = 3), 100)

eseq <- c(1, 2, 4, 8, 10, 40, 100)

eseq <- c(1, 2, 3, 4, 6, 8, 10, 16, 40, 100)


S <- 500
K <- 500
eseq <- c(7, 10, 20, 30, 70, 100, 200, 500)  # next: c(40, 60, 200)
qseq <- c(.05, .2)
qseq <- .05
cseq <- c(50, 200)
cseq <- 200

nseq  <- seq(0, 1.5*cc, by=4)
nsurv <- seq( 0, 1, length = 26 )
survMat <- matrix( 0, length(nsurv), length(nseq), dimnames = list(nsurv, nseq))

eseq <- rev(eseq)                        # must be descending

dseq <- c('beta','cov', 'unstr')

dseq <- 'unstr'

dseq <- 'beta'

dseq <- c('cov')


design <- expand.grid(eseq, dseq, qseq, cseq) 
colnames(design) <- c('E','DESIGN', 'q', 'C')
design$DESIGN <- as.character(design$DESIGN)
design$CONDDISP <- T
design$K <- K
design$S <- S

# ww <- which( design$DESIGN == 'unstr' )[1]
# design <- design[ 1:ww, ]


div <- matrix(0, nsim, 4)
colnames(div) <- c('diversity','richness','siteSorting','sppSorting') # landscape diversity
nexp <- nrow(design)

effects <- matrix(NA, nexp, 12)
ci <- c( 'Mu', paste( c('-','+'), 'SD', sep='' ) )
colnames(effects) <- as.vector( t(outer( colnames(div), ci, paste, sep='' )) )

survList <- vector( nexp, mode = 'list' )
names(survList) <- paste( design[,1], design[,2], sep = '_')
corNSMu <- numeric(0)


nt <- 300
nsim <- 5

ie <- 1
nm <- nexp

for(m in ie:nm){           # number of experiments
  
  nave  <- ntot <- matrix(0, nt, S)            # average over simulations
  betam <- beta
  DES   <- design$DESIGN[m]
  CONDDISP <- design$CONDDISP[m]
  
  if('DESIGN' %in% colnames(design)){
    DES <- design[m,'DESIGN']
  }
  
 # print( m )
  cat('\n\n')
  print( design[m,] )
  
  survDens <- survMat*0
  corNS <- numeric(0)
  
  E <- design[m,'E']
  q <- design[m,'q']
  cc <- design[m,'C']
  
  
  div <- matrix(0, nsim, 4)
  colnames(div) <- c('diversity','richness','siteSorting','sppSorting') # landscape diversity
  
  for(i in 1:nsim){  # repeat simulation
    
    beta <- matrix( rnorm(E*S, 0, sqrt(varBeta) ), E, S )     # X, B, V (matrix S in ms)
    vx   <- cov( rmvnormRcpp(10, rep(0, E), varX*diag(E) ) )
    x    <- rmvnormRcpp(K, rep(0, E), vx )
  #  x    <- sweep(x, 2, colMeans(x), '-')
    vx   <- cov(x)
    cs   <- t(beta)%*%vx%*%beta                                 # C = B'VB
    xb   <- x%*%beta                                            # XB
    colnames(x) <- rownames(beta) <- paste('q',1:nrow(beta), sep='')
    
    betam <- beta
    
    if( 'E' %in% colnames(design) ){
      q2    <- as.numeric( design[m,'E'] )
      betam <- beta[drop=F, 1:q2,]
      xb    <- x[,1:q2, drop=F]%*%betam
      vx    <- cov(x[,1:q2, drop=F])
      cs    <- t(betam)%*%vx%*%betam  # matrix C
    }
    
    
    c_k <- rowMeans(xb)
    c_s <- colMeans(xb)
    v_s <- diag( var(xb) )
    v_k <- diag( var(t(xb)) )
    
    competitiveSites <- cbind( c_k, sqrt(v_k) )
    competitiveSpp   <- cbind( c_s, sqrt(v_s) )
    colnames(competitiveSites) <- colnames(competitiveSpp) <- c('mean','sd')
    
  #  cat('\niteration: ')
    print(i)
    
    tmp <- simLoop( DES = DES, ninit, K, S, nt, cc, z, beta, cs, ss, tau,
                        q, xb, intercept = intercept, survMat = survMat, tstart = .75*nt, 
                    CONDDISP = CONDDISP )
    nmean <- tmp$nmean  # site by species
    ntot  <- tmp$ntot   # time by species
    survByDens <- tmp$survByDens
    
    byTot  <- colSums( round(tmp$nlast) )     # at least one individual
    rich   <- length( which(byTot > 1) )
    
  #  survDens[ rownames(survByDens), colnames(survByDens) ] <- 
  #    survDens[ rownames(survByDens), colnames(survByDens) ] + survByDens
    
    siteN <- rowMeans(nmean)
    siteS <- colMeans(nmean)
    corN <- cor(siteN, competitiveSites[,1])
    corS <- cor(siteS, competitiveSpp[,1])
    
    corNS <- rbind( corNS, c(corN, corS) )
    
    corNmean   <- suppressWarnings( cor(nmean, use = "complete.obs") )
    wf <- which( is.finite(corNmean) & !corNmean == 1 )
    crr <- cov2cor(cs) 
    sppSorting <- cor( as.vector( crr[wf] ), as.vector( corNmean[wf] ))
    
  #  plotTime( ntot )
    
    # sorting: site suitability vs mean abundance
    siteSorting <- cor( as.vector(xb), as.vector(nmean) )
    
    # species diversity
    tmp        <- diversity( nmean )
    div[i,] <- c( exp(tmp$byTot), rich, siteSorting, sppSorting )
    
    print(div[i,1:2])
    
    nave <- nave + ntot
    
  }       # replicate simlations
  
#  survDens 
  
 # sr <- range( which( rowSums(survDens) > 0 ) )
 # sc <- range( which( colSums(survDens) > 0 ) )
  
 # survByDens <- survDens[ sr[1]:sr[2], sc[1]:sc[2] ]
  
 # survList[[m]] <- survByDens
  corNSMu <- rbind( corNSMu, colMeans(corNS) )
  
  nave <- nave/nsim     # by-time average over simulations
  
  dnow <-  apply( div, 2, quantile, pnorm(c(0, -1, 1) ), na.rm=T )
  dnow[ is.na(dnow) ] <- 0
  
  effects[m,] <- as.vector( signif( dnow, 4 ) )
  
  print(effects[1:m,])
}  # design loop

effects <- data.frame(design, effects)
effects$E <- as.numeric(effects$E)

tmp <- read.csv("effects300x5.csv", stringsAsFactors = F)

tmp[,c(1:8, 11)]

############################# merge for 300x5

newColumns <- function(data){
  
  # add sd column
  
  kcols <- c("E","K", "q", "nsim", "DESIGN","CONDDISP")
  
  if( !'K' %in% colnames(data) )data$K <- 50
  if( !'nsim' %in% colnames(data) )dnew$nsim <- 5
  
  sc <- which( endsWith( colnames(data), '.SD') )
  sdcol <- data[,sc-1] - data[,sc]
  colnames(sdcol) <- .replaceString( colnames(sdcol), 'Mu', 'SD')
  
  dnew <- data[,kcols]
  
  dnew <- cbind( dnew, data[,sc-1], sdcol )
  ww <- grep('.SD', colnames(dnew), fixed=T )
  if(length(ww) > 0)dnew <- dnew[,-ww]
  dnew[,unique(colnames(dnew))]
}


enew <- newColumns( effects )
tnew <- newColumns( tmp )

anew <- rbind(enew, tnew[,colnames(enew) ] )

index <- paste( anew$DESIGN, anew$E, sep='-' )
nrows <- tapply( anew$nsim, index, sum )
mcols <- rep( grep('Mu',colnames(anew)), each = nrow(anew) )
arows <- rep( index, 4 )
mrows <- tapply( as.vector( unlist(anew[,grep('Mu',colnames(anew))]) ) , list( arows, mcols), mean )
colnames(mrows) <- colnames( anew )[grep('Mu',colnames(anew))]

srows <- tapply( as.vector( unlist(anew[,grep('SD',colnames(anew))]^2) ) , list( arows, mcols), sum )
srows <- sqrt( srows/2 )
colnames(srows) <- colnames( anew )[grep('SD',colnames(anew))]

irow <- sort( unique( index ) )
mi   <- match( irow, index )

atot <- anew[mi, c("E","DESIGN","q","C","CONDDISP","K","S")]
atot$nsim <- nrows[irow]
atot <- cbind( atot, mrows[irow,], srows[irow,] )
lo   <- mrows - 1.96*srows
hi   <- mrows + 1.96*srows
colnames(lo) <- .replaceString( colnames(lo) , 'Mu','Lo' )
colnames(hi) <- .replaceString( colnames(hi) , 'Mu','Hi' )
atot <- cbind( atot, lo, hi )
atot <- atot[ order( atot$DESIGN, atot$E), ]

write.csv( atot, file = 'effectsTot.csv', row.names=F )


e1 <- read.csv("effectsConddisp.csv", stringsAsFactors = F)
e1 <- e1[ e1$DESIGN %in% c('beta','cov'), ]
e1$nsim <- 10
e1$K <- 50
e1$q <- .1
e1$C <- 200
e1   <- e1[ e1$E > 1, ]
e1 <- newColumns( e1 )
mcol <- grep('Mu', colnames(e1))
scol <- which( endsWith(colnames(e1), 'SD') )
lo   <- e1[,mcol] - 1.96*e1[,scol]
hi   <- e1[,mcol] + 1.96*e1[,scol]

colnames(lo) <- .replaceString( colnames(lo) , 'Mu','Lo' )
colnames(hi) <- .replaceString( colnames(hi) , 'Mu','Hi' )
e1 <- cbind(e1, lo, hi )

effectsAll <- rbind( e1, atot[,colnames(e1)] )

# diversity


cramp <- colorRampPalette(mcolor)
cols <- cramp(3)

col2 <- colorRampPalette(qcolor)(3)

par(bty='n', mfrow=c(1,4), bty='n', mar=c(4,2,1,2), omi = c(.2,.1,.1,.1), cex = 1.01)

dleg <- c('SK', 'SU', 'UU')
vars <- c('diversity','richness', 'siteSorting', 'sppSorting')
ylabs <-  c('Diversity','Richness', 'Site Sorting', 'Species Sorting')
xlab <- ''
xlim <- c(1, 1100)
d    <- unique(atot$DESIGN)

for(k in 1:4){
  
  xlab <- ''
  mmm  <- paste(vars[k], 'Mu',sep='')
  mlo  <- paste(vars[k], 'Lo',sep='')
  mhi  <- paste(vars[k], 'Hi',sep='')
#  sdds <- paste(vars[k], c('.SD','.SD.1'), sep='')
  if( k %in% c(1,2) ){
    ylim <- c(5,200) 
    plot(NA, xlim = xlim, ylim = ylim, log = 'xy',
       xlab = xlab, ylab = ylabs[k])
  }else{
    ylim <- c(0, .7)
    plot(NA, xlim = xlim, ylim = ylim, log = 'x',
         xlab = xlab, ylab = ylabs[k])
  }
  if( k == 4 )abline(h=0, lty=2, lwd=2, col = 'grey')

  
  ni <- 3
  if(k == 4)ni <- 2
  
  for(i in 1:ni){
    
    wi <- which(atot$DESIGN == d[i])
    if( length(wi) > 0 ){
      yy   <- smooth( atot[wi,mmm], kind = '3' )
      lo   <- smooth( atot[wi, mlo], kind = '3' )
      hi   <- smooth( atot[wi, mhi], kind = '3' )
      lines( atot$E[wi], yy, col = cols[i], lwd=2)
      .shadeInterval( xvalues= atot$E[wi],
                      loHi = cbind( lo, hi ),
                      col= cols[i], add  = TRUE, trans = .4)
    }
    
    wi <- which(e1$DESIGN == d[i])
    yy   <- smooth(  e1[wi,mmm], kind = '3' )
    lo   <- smooth( e1[wi, mlo], kind = '3' )
    hi   <- smooth( e1[wi, mhi], kind = '3' )
    
    lines( e1$E[wi], yy, col = cols[i], lwd=2, lty=2)
    .shadeInterval( xvalues= e1$E[wi],
                    loHi = cbind( lo, hi ),
                    col= cols[i], add  = TRUE, trans = .4)
    
  #  if( k == 2 ){ 
    ushift <- -1.8
    if( k == 3 )ushift <- .075
      mi <- max(wi)
      wi <- max( which(atot$DESIGN == d[i]) )
      if( d[i] == 'beta' )text( 1.*e1$E[mi], e1[mi,mmm], 'SK', col = cols[i], pos=4 )
      if( d[i] == 'cov' )text( 1.*e1$E[mi], e1[mi,mmm], 'SU', col = cols[i], pos=4 )
      if( d[i] == 'unstr' )text( 1.2*atot$E[wi], ushift + atot[wi,mmm], 'UU', col = cols[i], pos=1 )
 #   }
  }
}

mtext('E', 1, outer=T, line = -1 )


##############################



tmp <- read.csv("effectsConddisp.csv", stringsAsFactors = F)
tmp <- read.csv("effectsNotConddisp.csv", stringsAsFactors = F)


enew <- rbind(effects, tmp[,colnames(effects)])
enew <- enew[ order(enew$CONDDISP, enew$DESIGN, enew$K, enew$S, enew$C, enew$q, enew$E ), ]
rownames(enew) <- NULL
enew[,c(1:8, 11)]

effects <- enew
write.csv( effects, file = "effects.csv", row.names = F)
write.csv( effects, file = "effects200x3.csv", row.names = F)



write.csv( enew, file = "effects300x5.csv", row.names = F)





greens <- c('#ffffcc','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32')
dcolor <- c('#d53e4f','#d53e4f','#fc8d59','#fee08b','#ffffbf','#e6f598','#99d594','#3288bd','#3288bd')
qcolor <- c('#7fc97f','#beaed4','#fdc086','#386cb0')
scolor <- c('white','#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#8c2d04')
bcolor <- c('white','#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84')
mcolor <- c('#66c2a5','#fc8d62','#8da0cb')
ocolor <- c('#ef6548','#7f0000')


################## surv 
par(mfrow=c(3,3), bty='n', mar = c(4,4,2,1))

for(k in 1:9){
  survByDens <- survList[[k]]
  survByDens <- survByDens[-1,]
  
  survProb <- sweep( survByDens, 2, colSums(survByDens), '/' ) # rows are surv, cols are abundance
  snames <- as.numeric( rownames(survProb) )
  nnames <- as.numeric( colnames(survProb) )
  survProb[ !is.finite(survProb) ] <- 0
  
  sn <- rep( snames, each = length(nnames))
  nn <- rep( nnames, length(snames))
  
  pseq   <- seq(0, max(survProb, na.rm=T), length=50)
  
  cramp <- colorRampPalette(bcolor)
  cols <- cramp(50)[ findInterval( (t(survProb))^.2, pseq, all.inside = T ) ]
  
  dx <- diff(nnames)[1] + nn*0
  dy <- diff(snames)[1] + sn*0
  
  symbols( x =nn, y = sn, rectangles = cbind(dx,dy), 
           inches = F, fg = cols, bg = cols,
           xlab = 'Local abundance', ylab = 'Survival probability', xlim = c(0, 200),
           ylim = c(0,1))
  title( paste0(design[k,], collapse = ' ' ) )
}
                                                  

write.csv( effects, file = 'effectsConddisp.csv', row.names = F )

row.names(effects) <- NULL



############## diversity



e1 <- read.csv("effectsConddisp.csv", stringsAsFactors = F)
e1 <- e1[ e1$DESIGN %in% c('beta','cov'), ]

e2 <- read.csv("effects300x5.csv", stringsAsFactors = F)
ww <- which( e2$CONDDISP & e2$C == 200 & e2$z == 100 & e2$q == .1 & e2$K == 50 & e2$S == 50 )
#e1 <- e1[ww, c(1:8, 11)]

sdc <- grep('.SD.1', colnames(e1))
muc <- grep('Mu', colnames(e1))
sd  <- e1[,sdc] - e1[,muc]
lo  <- e1[,muc] - 1.96*sd
hi  <- e1[,muc] + 1.96*sd
colnames(lo) <- sub('Mu','Lo',colnames(lo))
colnames(hi) <- sub('Mu','Hi',colnames(hi))
f1 <- cbind( e1[, c(1:2, muc )], lo, hi )

sdc <- grep('.SD.1', colnames(e2))
muc <- grep('Mu', colnames(e2))
sd  <- e2[,sdc] - e2[,muc]
lo  <- e2[,muc] - 1.96*sd
hi  <- e2[,muc] + 1.96*sd
colnames(lo) <- sub('Mu','Lo',colnames(lo))
colnames(hi) <- sub('Mu','Hi',colnames(hi))
f2 <- cbind( e2[, c(1:2, muc )], lo, hi )


cramp <- colorRampPalette(mcolor)
cols <- cramp(3)

col2 <- colorRampPalette(qcolor)(3)

par(bty='n', mfrow=c(1,4), bty='n', mar=c(4,3,1,2), omi = c(.2,.1,.1,.1), cex = 1.01)

dleg <- c('SK', 'SU', 'UU')
vars <- c('diversity','richness', 'siteSorting', 'sppSorting')
ylabs <-  c('Diversity','Richness', 'Site Sorting', 'Species Sorting')
xlab <- ''
xlim <- c(2, 1000)

for(k in 1:4){
  
  ylim <- c(0, .6)
  if( vars[k] %in% c('diversity','richness'))ylim <- c(0,150) 
  
  xlab <- ''
  mmm  <- paste(vars[k], 'Mu',sep='')
  mlo  <- paste(vars[k], 'Lo',sep='')
  mhi  <- paste(vars[k], 'Hi',sep='')
  sdds <- paste(vars[k], c('.SD','.SD.1'), sep='')
  plot(NA, xlim = xlim, ylim = ylim, log = 'x',
       xlab = xlab, ylab = ylabs[k])
  if( k == 4 )abline(h=0, lty=2, lwd=2, col = 'grey')
  d <- unique(e2$DESIGN)
  
  ni <- 3
  if(k == 4)ni <- 2
  for(i in 1:ni){
    
    wi <- which(e2$DESIGN == d[i])
    if(length(wi) > 0){
      lines( f2$E[wi], e2[wi,mmm], col = cols[i], lwd=2)
      .shadeInterval( xvalues= f2$E[wi],
                      loHi = f2[wi, c(mlo, mhi)],
                      col= cols[i], add  = TRUE, trans = .4)
    }
    
    wi <- which(e1$DESIGN == d[i])
    lines( e1$E[wi], e1[wi,mmm], col = cols[i], lwd=2, lty=2)
    .shadeInterval( xvalues= e1$E[wi],
                    loHi = e1[wi, sdds],
                    col= cols[i], add  = TRUE, trans = .4)
    
    if( k == 2 ){ 
      mi <- max(wi)
      if( d[i] == 'beta' )text( 1.1*e1$E[mi], e1[mi,mmm], 'SK', col = cols[i], pos=4 )
      if( d[i] == 'cov' )text( 1.1*e1$E[mi], e1[mi,mmm], 'SU', col = cols[i], pos=4 )
      if( d[i] == 'unstr' )text( 1.1*e1$E[mi], e1[mi,mmm], 'UU', col = cols[i], pos=4 )
    }
  }
}

#legend('bottomright', dleg, text.col = cols, bty='n' )
mtext('E', 1, outer=T, line = -1 )


write.csv(effects, file = 'effects.csv', row.names = F)


# compare CONDISP

cdisp <- read.csv( 'effectsConddisp.csv', stringsAsFactors = F )
cdisp <- cdisp[ cdisp$DESIGN == 'cov', ]
ndisp <- read.csv( 'effectsNotConddisp.csv', stringsAsFactors = F )
ndisp <- ndisp[ ndisp$DESIGN == 'cov', ]

cols <- ocolor

par(bty='n', mfrow=c(1,4), omi = c(.5,.3,.1,.1), bty='n')

for(k in 1:4){
  
  ylim <- range(effects[,grep(vars[k], colnames(effects))])
  if( vars[k] %in% c('diversity','richness')){
    ylim <- c(0,100) 
  }else if(vars[k] == 'siteSorting'){
    ylim <- c(-.1, .1)
  }else{
    ylim <- c(0, .7)
  }
  
  xlab <- ''
  mmm  <- paste(vars[k], 'Mu',sep='')
  sdds <- paste(vars[k], c('.SD','.SD.1'), sep='')
  
  nn <- 1:length(xc)
  if( k == 4 )nn <- nn[-1]
  
  plot(NA, xlim = range(effects$E), 
       ylim = ylim, 
       xlab = xlab, ylab = vars[k])
  lines( cdisp$E[nn], cdisp[nn,mmm], col = cols[1], lwd = 2 )
  .shadeInterval( cdisp$E[nn], cdisp[nn,sdds], col = cols[1], add=T, trans = .4)
  
  lines( ndisp$E[nn], ndisp[nn,mmm], col = cols[2], lwd = 2 )
  .shadeInterval( ndisp$E[nn], ndisp[nn,sdds], col = cols[2], add=T, trans = .4)
}
legend('topright', c('conditional dispersal','unconditional'), text.col = cols, bty='n' )
mtext('E', 1, outer=T )




title( DESIGN )

mtext('Time', side=1, line=.4, outer  = TRUE, cex=1.2)
mtext('Abundance', side=2, line=.4, outer  = TRUE, cex=1.2)
mtext( paste( 'E = ', design[1], sep=''), side=3, line=.4, outer  = TRUE, cex=1.2)
mtext('E = 2', side=3, line=-17, outer  = TRUE, cex=1.2)





# species covariance

par(bty='n', mfrow=c(1,2), mar=c(1,1,1,1))

cramp <- colorRampPalette( rev(dcolor) )
cramp <- cramp(40)

#cmm  <- max( abs(cs) )
#cseq <- seq(-cmm, cmm, length=40)

q2 <- 2

x     <- rmvnormRcpp(50, rep(0, E), diag(E) )
beta  <- matrix( rnorm(E*S, 0, sqrt(varBeta) ), E, S )
betam <- beta[1:q2,]
xb    <- x[,1:q2]%*%betam
vx    <- cov(x[,1:q2])
cs    <- t(betam)%*%vx%*%betam  # matrix C
#cs    <- cov2cor(cs)


cmm  <- max( abs(cs) )
cseq <- seq(-cmm, cmm, length=40)

cols <- findInterval( cs, cseq, all.inside = T)
cols <- cramp[ cols ]
ix   <- expand.grid( 1:S, 1:S )
wx   <- which(ix[,1] >= ix[,2])

plot( ix[wx,1], ix[wx,2], col = cols[wx], xaxt='n', yaxt='n',
      xlab = '', ylab = '', pch=15, cex=1, asp=1)
ch   <- cs
ch[ which(lower.tri(ch)) ] <- -99
kmax <- apply(ch, 2, which.max )
ik   <- cbind( 1:S, kmax )
wx <- which(ik[,2] <= ik[,1])
#points(ik[,2], ik[,1], pch=0, cex=.5, col='black')
points(ik[wx,1], ik[wx,2], pch=0, cex=1, col='black')

title( paste('E =', q2, ', S =', S  ))

q2 <- 20
betam <- beta[1:q2,]
xb    <- x[,1:q2]%*%betam
vx    <- cov(x[,1:q2])
cs    <- t(betam)%*%vx%*%betam  # matrix C
#cs    <- cov2cor(cs)
cmm  <- max( abs(cs) )
cseq <- seq(-cmm, cmm, length=40)

print(cmm)

cols <- findInterval( cs, cseq, all.inside = T)
cols <- cramp[ cols ]
ix   <- expand.grid( 1:S, 1:S )
wx   <- which(ix[,1] >= ix[,2])

plot( ix[wx,1], ix[wx,2], col = cols[wx], xaxt='n', yaxt='n',
      xlab = '', ylab = '', pch=15, cex=1, asp=1)
ch   <- cs
ch[ which(lower.tri(ch)) ] <- -99
kmax <- apply(ch, 2, which.max )
ik   <- cbind( 1:S, kmax )
wx <- which(ik[,2] <= ik[,1])
#points(ik[,2], ik[,1], pch=0, cex=.5, col='black')
points(ik[wx,1], ik[wx,2], pch=0, cex=1, col='black')

title( paste('E = ', q2, ', S =', S  ))








# outcomes map
#cc  <- 5*10        
#z  <- cc/2 
q2 <- 50
S  <- 100
E  <- 500

vx   <- cov( rmvnormRcpp(10, rep(0, E), varX*diag(E) ) )
x    <- rmvnormRcpp(K, rep(0, E), vx )
vx   <- cov(x)
xb    <- x[,1:q2]%*%betam
beta  <- matrix( rnorm(E*S, 0, sqrt(varBeta) ), E, S )
betam <- beta[drop=F, 1:q2,]
xb    <- x[,1:q2, drop=F]%*%betam
vx    <- cov(x[,1:q2, drop=F])
cs    <- t(betam)%*%vx%*%betam  # matrix C


tmp <- simLoop(DES = 'beta', ninit, K, S, nt, cc, z, beta = betam, cs, ss, tau,
               q, xb, intercept = intercept, survMat = survMat, tstart = .75*nt, 
               CONDDISP = T )

tmp <- simLoop(DES = 'cov', ninit, K, S, nt, cc, z, beta = betam, cs, ss, tau,
               q, xb, intercept = intercept, survMat = survMat, tstart = .75*nt, 
               CONDDISP = T )

nlast <- tmp$nlast
  
cramp <- colorRampPalette(greens)
cramp <- cramp(40)
cseq <- seq(min(xb), max(xb), length=40)
cols <- findInterval( xb, cseq, all.inside = T)
cols <- cramp[ cols ]
ix   <- expand.grid( 1:S, 1:K )

par(bty='n', mfrow=c(1,2))
plot( ix[,1], ix[,2], col = cols, xaxt='n', yaxt='n',
      xlab = 'Species', ylab = '', pch=15, cex=.6, asp=1)
kmax <- apply(xb, 1, which.max )
ik   <- cbind( 1:K, kmax )
points(ik[,1], ik[,2], pch=0, cex=.6)

plot( xb, nlast, xlab='Suitability', ylab='', pch = 0)#, log = 'y' )
points( xb[ ik ], nlast[ ik ], col = 'darkgreen', pch=15)

mtext('Site', side=2, line=-5, outer  = TRUE, cex=1.2)



