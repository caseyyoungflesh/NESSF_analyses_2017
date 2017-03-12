######################
#Characterizing change in prediced dN over time at Adelie penguin colonies
#
#Script to convolve lab collected guano spectra to Landsat bands
#
#Authors: Casey Youngflesh
######################



# Load guano spectra and SIA results -----------------------------------

#spectra from guano samples collected in the field
setwd(paste0(dir, 'Master_spectra'))

t2 <- read.csv('Wet_processed_spectra_14_15.csv')


#merge with SIA data
#setwd('/Users/caseyyoungflesh/Google Drive/R/Guano_SIA/')
setwd('C:/Users/Lynch Lab 7/Google Drive/R/Guano_SIA/')

SIA_data <- read.csv('SIA_lsu.csv', header=TRUE)
names(SIA_data)[4] <- 'ID'

#wet spec merged with SIA data
w_data <- merge(SIA_data, t2, by='ID')


# Function to convolve spectra to sensor bands ----------------------------


#Select sensor L7, L8, WV2, WV3
#Data is WL in columns, samples in rows
#Data can have other columns before WL start
#Starting column for spec info is argument 'stcol'


convolve_spectra <- function(Data, sensor= 'L7', stcol=24)
{
  #Data <- w_data
  #subset to not include info before WL dat
  DAT <- Data[,stcol:NCOL(Data)]
  
  #####
  #Begin bands
  #####
  #Landsat 4/5 bands
  L45.band1 <- c(450, 520)
  L45.band2 <- c(520, 600)
  L45.band3 <- c(630, 690)
  L45.band4 <- c(760, 900)
  L45.band5 <- c(1550, 1750)
  #L45.band6 <- c(10400, 12500) #thermal band
  L45.band7 <- c(2080, 2350)
  
  
  #Landsat 7 bands
  L7.band1 <- c(450, 520)
  L7.band2 <- c(520, 600)
  L7.band3 <- c(630, 690)
  L7.band4 <- c(770, 900)
  L7.band5 <- c(1550, 1750)
  #L7.band6 <- c(10400, 12500) #thermal band
  L7.band7 <- c(2090, 2350)
  #L7.band8 <- c(520, 900) #panchromatic
  
  #Landsat 8 bands
  L8.band1 <- c(430, 450)
  L8.band2 <- c(450, 510)
  L8.band3 <- c(530, 590)
  L8.band4 <- c(640, 670)
  L8.band5 <- c(850, 880)
  L8.band6 <- c(1570, 1650)
  L8.band7 <- c(2110, 2290)
  #L8.band8 <- c(500, 680) #panchromatic
  #L8.band9 <- c(1360, 1380) #cirrus
  #L8.band10 <- c(10600, 11190) #thermal1
  #L8.band11 <- c(11500, 12510) #thermal2
  
  #WordlsView 2 bands
  WV2.band1 <- c(400, 450)
  WV2.band2 <- c(450, 510) 
  WV2.band3 <- c(510, 580)
  WV2.band4 <- c(585, 625) 
  WV2.band5 <- c(630, 690)
  WV2.band6 <- c(705, 745) 
  WV2.band7 <- c(770, 895)
  WV2.band8 <- c(860, 1040) 
  #WV2.band9 <- c(450, 800) #panchromatic
  
  #WordlsView 3 bands
  WV3.band1 <- c(400, 450)
  WV3.band2 <- c(450, 510) 
  WV3.band3 <- c(510, 580)
  WV3.band4 <- c(585, 625) 
  WV3.band5 <- c(630, 690)
  WV3.band6 <- c(705, 745) 
  WV3.band7 <- c(770, 895)
  WV3.band8 <- c(860, 1040) 
  #WV3.band9 <- c(450, 800) #panchromatic
  WV3.band10 <- c(1195, 1225)
  WV3.band11 <- c(1550, 1590) 
  WV3.band12 <- c(1640, 1680)
  WV3.band13 <- c(1710, 1750) 
  WV3.band14 <- c(2145, 2185)
  WV3.band15 <- c(2185, 2225) 
  WV3.band16 <- c(2235, 2285)
  WV3.band17 <- c(2295, 2365)
  #####
  #End bands
  #####
  
  #put other data info in object OTHER
  if(stcol > 1)
  {
    OTHER <- Data[,1:stcol-1]
  }
  
  
  if(sensor == 'L45')
  {
    temp2<- c()
    for (i in 1:NROW(DAT))
    {
      #i <- 1
      
      OUT1 <- mean(as.numeric(DAT[i,(L45.band1-350)[1]:(L45.band1-350)[2]]))
      OUT2 <- mean(as.numeric(DAT[i,(L45.band2-350)[1]:(L45.band2-350)[2]]))
      OUT3 <- mean(as.numeric(DAT[i,(L45.band3-350)[1]:(L45.band3-350)[2]]))
      OUT4 <- mean(as.numeric(DAT[i,(L45.band4-350)[1]:(L45.band4-350)[2]]))
      OUT5 <- mean(as.numeric(DAT[i,(L45.band5-350)[1]:(L45.band5-350)[2]]))
      OUT7 <- mean(as.numeric(DAT[i,(L45.band7-350)[1]:(L45.band7-350)[2]]))
      temp <- cbind(OUT1, OUT2, OUT3, OUT4, OUT5, OUT7)
      temp2 <- rbind(temp2, temp)
    }
    
    if(stcol > 1)
    {
      colnames(temp2) <- c('L45.Band1', 'L45.Band2', 'L45.Band3', 'L45.Band4', 'L45.Band5', 'L45.Band7')
      OUT.final <- cbind(OTHER, temp2)
    }else{
      OUT.final <- temp2
    }
  }
  
  
  
  
  if(sensor == 'L7')
  {
    temp2<- c()
    for (i in 1:NROW(DAT))
    {
      #i <- 1
      
      OUT1 <- mean(as.numeric(DAT[i,(L7.band1-350)[1]:(L7.band1-350)[2]]))
      OUT2 <- mean(as.numeric(DAT[i,(L7.band2-350)[1]:(L7.band2-350)[2]]))
      OUT3 <- mean(as.numeric(DAT[i,(L7.band3-350)[1]:(L7.band3-350)[2]]))
      OUT4 <- mean(as.numeric(DAT[i,(L7.band4-350)[1]:(L7.band4-350)[2]]))
      OUT5 <- mean(as.numeric(DAT[i,(L7.band5-350)[1]:(L7.band5-350)[2]]))
      OUT7 <- mean(as.numeric(DAT[i,(L7.band7-350)[1]:(L7.band7-350)[2]]))
      temp <- cbind(OUT1, OUT2, OUT3, OUT4, OUT5, OUT7)
      temp2 <- rbind(temp2, temp)
    }
    
    if(stcol > 1)
    {
      colnames(temp2) <- c('L7.Band1', 'L7.Band2', 'L7.Band3', 'L7.Band4', 'L7.Band5', 'L7.Band7')
      OUT.final <- cbind(OTHER, temp2)
    }else{
      OUT.final <- temp2
    }
  }
  
  if(sensor == 'L8')
  {
    temp2<- c()
    for (i in 1:NROW(DAT))
    {
      #i <- 1
      
      OUT1 <- mean(as.numeric(DAT[i,(L8.band1-350)[1]:(L8.band1-350)[2]]))
      OUT2 <- mean(as.numeric(DAT[i,(L8.band2-350)[1]:(L8.band2-350)[2]]))
      OUT3 <- mean(as.numeric(DAT[i,(L8.band3-350)[1]:(L8.band3-350)[2]]))
      OUT4 <- mean(as.numeric(DAT[i,(L8.band4-350)[1]:(L8.band4-350)[2]]))
      OUT5 <- mean(as.numeric(DAT[i,(L8.band5-350)[1]:(L8.band5-350)[2]]))
      OUT6 <- mean(as.numeric(DAT[i,(L8.band6-350)[1]:(L8.band6-350)[2]]))
      OUT7 <- mean(as.numeric(DAT[i,(L8.band7-350)[1]:(L8.band7-350)[2]]))
      temp <- cbind(OUT1, OUT2, OUT3, OUT4, OUT5, OUT6, OUT7)
      temp2 <- rbind(temp2, temp)
    }
    
    if(stcol > 1)
    {
      colnames(temp2) <- c('L8.Band1', 'L8.Band2', 'L8.Band3', 'L8.Band4', 'L8.Band5', 
                           'L8.Band6', 'L8.Band7')
      OUT.final <- cbind(OTHER, temp2)
    }else{
      OUT.final <- temp2
    }
  }
  
  if(sensor == 'WV2')
  {
    temp2<- c()
    for (i in 1:NROW(DAT))
    {
      #i <- 1
      
      OUT1 <- mean(as.numeric(DAT[i,(WV2.band1-350)[1]:(WV2.band1-350)[2]]))
      OUT2 <- mean(as.numeric(DAT[i,(WV2.band2-350)[1]:(WV2.band2-350)[2]]))
      OUT3 <- mean(as.numeric(DAT[i,(WV2.band3-350)[1]:(WV2.band3-350)[2]]))
      OUT4 <- mean(as.numeric(DAT[i,(WV2.band4-350)[1]:(WV2.band4-350)[2]]))
      OUT5 <- mean(as.numeric(DAT[i,(WV2.band5-350)[1]:(WV2.band5-350)[2]]))
      OUT6 <- mean(as.numeric(DAT[i,(WV2.band6-350)[1]:(WV2.band6-350)[2]]))
      OUT7 <- mean(as.numeric(DAT[i,(WV2.band7-350)[1]:(WV2.band7-350)[2]]))
      OUT8 <- mean(as.numeric(DAT[i,(WV2.band8-350)[1]:(WV2.band8-350)[2]]))
      temp <- cbind(OUT1, OUT2, OUT3, OUT4, OUT5, OUT6, OUT7, OUT8)
      temp2 <- rbind(temp2, temp)
    }
    
    if(stcol > 1)
    {
      colnames(temp2) <- c('WV2.Band1', 'WV2.Band2', 'WV2.Band3', 'WV2.Band4', 'WV2.Band5', 
                           'WV2.Band6', 'WV2.Band7', 'WV2.Band8')
      OUT.final <- cbind(OTHER, temp2)
    }else{
      OUT.final <- temp2
    }
  }
  
  if(sensor == 'WV3')
  {
    temp2<- c()
    for (i in 1:NROW(DAT))
    {
      #i <- 1
      
      OUT1 <- mean(as.numeric(DAT[i,(WV3.band1-350)[1]:(WV3.band1-350)[2]]))
      OUT2 <- mean(as.numeric(DAT[i,(WV3.band2-350)[1]:(WV3.band2-350)[2]]))
      OUT3 <- mean(as.numeric(DAT[i,(WV3.band3-350)[1]:(WV3.band3-350)[2]]))
      OUT4 <- mean(as.numeric(DAT[i,(WV3.band4-350)[1]:(WV3.band4-350)[2]]))
      OUT5 <- mean(as.numeric(DAT[i,(WV3.band5-350)[1]:(WV3.band5-350)[2]]))
      OUT6 <- mean(as.numeric(DAT[i,(WV3.band6-350)[1]:(WV3.band6-350)[2]]))
      OUT7 <- mean(as.numeric(DAT[i,(WV3.band7-350)[1]:(WV3.band7-350)[2]]))
      OUT8 <- mean(as.numeric(DAT[i,(WV3.band8-350)[1]:(WV3.band8-350)[2]]))
      OUT10 <- mean(as.numeric(DAT[i,(WV3.band10-350)[1]:(WV3.band10-350)[2]]))
      OUT11 <- mean(as.numeric(DAT[i,(WV3.band11-350)[1]:(WV3.band11-350)[2]]))
      OUT12 <- mean(as.numeric(DAT[i,(WV3.band12-350)[1]:(WV3.band12-350)[2]]))
      OUT13 <- mean(as.numeric(DAT[i,(WV3.band13-350)[1]:(WV3.band13-350)[2]]))
      OUT14 <- mean(as.numeric(DAT[i,(WV3.band14-350)[1]:(WV3.band14-350)[2]]))
      OUT15 <- mean(as.numeric(DAT[i,(WV3.band15-350)[1]:(WV3.band15-350)[2]]))
      OUT16 <- mean(as.numeric(DAT[i,(WV3.band16-350)[1]:(WV3.band16-350)[2]]))
      OUT17 <- mean(as.numeric(DAT[i,(WV3.band17-350)[1]:(WV3.band17-350)[2]]))
      
      
      temp <- cbind(OUT1, OUT2, OUT3, OUT4, OUT5, OUT6, OUT7, OUT8, OUT10,
                    OUT11, OUT12, OUT13, OUT14, OUT15, OUT16, OUT17)
      temp2 <- rbind(temp2, temp)
    }
    
    if(stcol > 1)
    {
      colnames(temp2) <- c('WV3.Band1', 'WV3.Band2', 'WV3.Band3', 'WV3.Band4', 'WV3.Band5', 
                           'WV3.Band6', 'WV3.Band7', 'WV3.Band8', 'WV3.Band10', 'WV3.Band11'
                           , 'WV3.Band12', 'WV3.Band13', 'WV3.Band14', 'WV3.Band15', 'WV3.Band16'
                           , 'WV3.Band17')
      OUT.final <- cbind(OTHER, temp2)
    }else{
      OUT.final <- temp2
    }
  }
  
  return(OUT.final)
}







# Convolve L45, L7, and L8 spectra using function -------------------------

conv_L45 <- convolve_spectra(w_data, sensor = 'L45', stcol = 24)
conv_L7 <- convolve_spectra(w_data, sensor = 'L7', stcol = 24)
conv_L8 <- convolve_spectra(w_data, sensor = 'L8', stcol = 24)


