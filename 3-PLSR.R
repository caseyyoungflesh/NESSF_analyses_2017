######################
#Characterizing change in prediced dN over time at Adelie penguin colonies
#
#Script to examine relationship between convolve lab spectra and dN as obtained by SIA
#
#Authors: Casey Youngflesh
######################




# Functions for PLSR analysis ----------------------------------------------

#-------------#
#Function to determine optimal number of components

#colony data (response) and colony covariates (predictors)
#comps = number of components
#iterations = number of iterations

comp.fun <- function(COL_DATA, COL_COV, comps=10, iterations=100)
{
  
  #COL_COV <- bands
  #COL_DATA <- J_YEAR
  
  # Jackknife test here.  Determine optimal number of components
  dims = dim(COL_COV)
  
  jk.out <- array(data=NA,dim=c(iterations,comps+1)) # num of components plus intercept
  for (i in 1:iterations)
  {
    #i <- 1
    print(paste("Iteration: ",i,sep=""))
    #segs = cvsegments(dims[1],k = k, type="random")
    trun.pls = plsr(COL_DATA ~ COL_COV, scale=FALSE, ncomp=comps, validation="CV",
                    trace=TRUE, method = "oscorespls")
    vec <- as.vector(RMSEP(trun.pls)$val)
    vec.sub.adj <- vec[seq(2,length(vec),2)]
    jk.out[i,] <- vec.sub.adj
  }
  
  temp <- deparse(substitute(COL_DATA))
  
  # Boxplot of results
  par(mfrow=c(1,1))
  numcomps <- comps+1
  boxplot(jk.out[,1:numcomps], xaxt="n", xlab="NUM OF COMPONENTS", ylab="RMSEP",
          main= paste(temp))
  axis(1,at=1:numcomps, 0:comps)
  box(lwd=2.2)
  
  return(jk.out)
}


#-------------#
#function to look at optimal # comps and see how well plsr performs

vis.fun <- function(penguin, cov)
{
  #penguin <- COPA.data$LAY
  #cov <- COPA.cov
  temp <- plsr(penguin ~ cov, validation='CV')
  par(mfrow=c(1,2))
  comp.fun(penguin, cov, comps=10)
  plot(R2(temp))
  
  return(summary(temp))
}

#-------------#
#VIP function

VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


#-------------#
#Function ONE to run

#determes # of components and variance explained
ONE <- function(RESPONSE, PREDICTOR, COMP = 10, IT = 100)
{
  #RESPONSE <- unlist(m_data$d_15N)
  #PREDICTOR <- samp.wv
  temp <- plsr(RESPONSE ~ PREDICTOR, validation='CV', method= 'oscorespls')
  par(mfrow=c(1,1))
  comp.fun(RESPONSE, PREDICTOR, comps= COMP, iterations = IT)
  return(summary(temp))
}


#-------------#
#Function TWO to run

#plots prediction plot (prediction vs measured) COEF estimate and VIP
TWO <- function(RESPONSE, PREDICTOR, COMPS)
{
  .pardefault <- graphics::par(no.readonly = T)
  tq <- plsr(RESPONSE ~ PREDICTOR, validation='CV', method= 'oscorespls')
  
  predplot(tq, ncomp = COMPS, pch=19, asp = 1, line = TRUE, main = 'Prediction Plot')
  
  par(mfrow=c(2,1))
  #COEF
  coefs = coef(tq, ncomp= COMPS,intercept=FALSE) # WITHOUT INTERCEPT FOR PLOTTING
  plot(coefs, xlab="Wavelength (nm)",ylab="REG COEF", pch=20, xaxt= 'n', main='COEF')
  axis(1, at= seq(151,2151, by=250), labels= seq(500,2500, by=250))
  lines(coefs,lwd=2.5)
  abline(h=0,lty=2,col="dark grey")
  
  #VIP
  vips = VIP(tq)[COMPS,]
  plot(vips,xlab="Wavelength (nm)",ylab="VIP", pch=20, xaxt = 'n', main='VIP')
  axis(1, at= seq(151,2151, by=250), labels= seq(500,2500, by=250))
  lines(vips,lwd=3)
  abline(h=0.8,lty=2,col="dark grey")
  graphics::par(.pardefault)
}





# PLSR Landsat-convolved lab spec to SIA ----------------------------------


#L45

#dC
response <- as.numeric(unlist(conv_L45$d_13C))
predictor <- data.matrix(conv_L45[,-c(1:23)])
ONE(response, predictor, COMP = 6, IT = 20) #4 comps = 63% variation

#dN
response <- as.numeric(unlist(conv_L45$d_15N))
predictor <- data.matrix(conv_L45[,-c(1:23)])
ONE(response, predictor, COMP = 6, IT = 20) #6 comps = 76% variation




#L7

#dC
response <- as.numeric(unlist(conv_L7$d_13C))
predictor <- data.matrix(conv_L7[,-c(1:23)])
ONE(response, predictor, COMP = 6, IT = 20) #4 comps = 63% variation

#dN
response <- as.numeric(unlist(conv_L7$d_15N))
predictor <- data.matrix(conv_L7[,-c(1:23)])
ONE(response, predictor, COMP = 6, IT = 20) #6 comps = 76% variation




#L8

#dC
response <- as.numeric(unlist(conv_L8$d_13C))
predictor <- data.matrix(conv_L8[,-c(1:23)])
ONE(response, predictor, COMP = 7, IT = 20) #7 comps - 73% variation

#dN
response <- as.numeric(unlist(conv_L8$d_15N))
predictor <- data.matrix(conv_L8[,-c(1:23)])
ONE(response, predictor, COMP = 7, IT = 20) #7 comps - 82% variation

