######################
#Characterizing change in prediced dN over time at Adelie penguin colonies
#
#Script to process data
#
#Authors: Casey Youngflesh
######################



# Clear environment -------------------------------------------------------

rm(list = ls())
dev.off()



# Set wd ------------------------------------------------------------------

osx <- '/Users/caseyyoungflesh/Google Drive/R/Antarctica_spec/'
win <- 'C:/Users/Lynch Lab 7/Google Drive/R/Antarctica_spec/'

if(Sys.info()[['sysname']] == 'Windows')
{
  dir <- win
}
if(Sys.info()[['sysname']] == 'Darwin')
{
  dir <- osx
}



# Load packages -----------------------------------------------------------

#installs packages if needed and loads required packages

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages(pacman)
}

pacman::p_load(dplyr, pls, R2jags, MCMCvis, ggplot2, mapproj)



# Load data ---------------------------------------------------------------

setwd(paste0(dir, 'Spec_change_over_time_Schwaller_algo/Data'))

L45 <- read.csv('L45_retrieval_dec_13_2016.csv', header = TRUE)
L7 <- read.csv('L7_retrieval_dec_13_2016.csv', header = TRUE)
L8 <- read.csv('L8_retrieval_dec_13_2016.csv', header = TRUE)

#Landsat file naming convention:
#LXSPPPRRRYYYYDDDGSIVV
#L = Landsat
#X = Sensor
#S = Satellite
#PPP = WRS path
#RRR = WRS row
#YYYY= year
#DDD = Julian day of year
#GSI = Ground station id
#VV = Archive version number



# Process L45 -------------------------------------------------------------

L45_OUT <- c()
for (i in 1:nrow(L45))
{
  #i <- 1
  YEAR <- as.numeric(substr(L45$source_fil[i], 10, 13))
  DAY <- as.numeric(substr(L45$source_fil[i], 14, 16))
  temp <- cbind(YEAR, DAY)
  L45_OUT <- rbind(L45_OUT, temp)
}


L45_1 <- mutate(L45, YEAR = L45_OUT[,1], 
                DAY = L45_OUT[,2], 
                SAT = rep('L45', nrow(L45)))


L45_p <- select(L45_1, SAT, CODE = X4_letter_c, 
                YEAR, DAY, D_VAL = d_value,
                SOURCE = source_fil,
                latitude, longitude, 
                coastal, blue, green, 
                red, NIR, SWIR1, SWIR2, cirrus)



# Process L7 --------------------------------------------------------------


L7_OUT <- c()
for (i in 1:nrow(L7))
{
  #i <- 1
  YEAR <- as.numeric(substr(L7$source_fil[i], 10, 13))
  DAY <- as.numeric(substr(L7$source_fil[i], 14, 16))
  temp <- cbind(YEAR, DAY)
  L7_OUT <- rbind(L7_OUT, temp)
}


L7_1 <- mutate(L7, YEAR = L7_OUT[,1], 
               DAY = L7_OUT[,2], 
               SAT = rep('L7', nrow(L7)))


L7_p <- select(L7_1, SAT, CODE = X4_letter_c, 
               YEAR, DAY, D_VAL = d_value,
               SOURCE = source_fil,
               latitude, longitude, 
               coastal, blue, green, 
               red, NIR, SWIR1, SWIR2, cirrus)



# Process L8 --------------------------------------------------------------


L8_OUT <- c()
for (i in 1:nrow(L8))
{
  #i <- 1
  YEAR <- as.numeric(substr(L8$source_fil[i], 10, 13))
  DAY <- as.numeric(substr(L8$source_fil[i], 14, 16))
  temp <- cbind(YEAR, DAY)
  L8_OUT <- rbind(L8_OUT, temp)
}
str(L8)

L8_1 <- mutate(L8, YEAR = L8_OUT[,1], 
               DAY = L8_OUT[,2], 
               SAT = rep('L8', nrow(L8)))


L8_p <- select(L8_1, SAT, CODE = X4_letter_c, 
               YEAR, DAY, D_VAL = d_value,
               SOURCE = source_fil,
               latitude, longitude, 
               coastal, blue, green, 
               red, NIR, SWIR1, SWIR2, cirrus)



# combine into single data.frame -----------------------------------------------------

data <- rbind(L45_p, L7_p, L8_p)

