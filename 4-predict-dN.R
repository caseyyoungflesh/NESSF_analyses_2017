######################
#Characterizing change in prediced dN over time at Adelie penguin colonies
#
#Script to predcit dN values based on PLSR analyses
#
#Authors: Casey Youngflesh
######################



# Predict dN from PLSR ------------------------

#predicts dN for each individual pixel separately

#Landsat 4/5
#L45 dN
response <- as.numeric(unlist(conv_L45$d_15N))
predictor <- data.matrix(conv_L45[,-c(1:23)])
L45_dN_fit <- plsr(response ~ predictor, ncomp = 6, validation='CV', method= 'oscorespls')


#predict new dN from Landsat 45

#just L45 retrievals - all pixels
L45_ret <- filter(data, SAT == 'L45')
L45_ret_bands <- select(L45_ret, blue, green, red, NIR, SWIR1, SWIR2)/10000

pred_fit <- predict(L45_dN_fit, ncomp = 6, newdata = data.matrix(L45_ret_bands))

#combine predictions with input data
L45_comb <- mutate(L45_ret, PRED= as.numeric(unlist(pred_fit)))



#Landsat 7
#L7 dN
response <- as.numeric(unlist(conv_L7$d_15N))
predictor <- data.matrix(conv_L7[,-c(1:23)])
L7_dN_fit <- plsr(response ~ predictor, ncomp = 6, validation='CV', method= 'oscorespls')


#predict new dN from Landsat 7

#just L7 retrievals - all pixels
L7_ret <- filter(data, SAT == 'L7')
L7_ret_bands <- select(L7_ret, blue, green, red, NIR, SWIR1, SWIR2)/10000

pred_fit <- predict(L7_dN_fit, ncomp = 6, newdata = data.matrix(L7_ret_bands))

#combine predictions with input data
L7_comb <- mutate(L7_ret, PRED= as.numeric(unlist(pred_fit)))



#Landsat 8
#L8 dN
response <- as.numeric(unlist(conv_L8$d_15N))
predictor <- data.matrix(conv_L8[,-c(1:23)])
L8_wet_fit <- plsr(response ~ predictor, ncomp = 7, validation='CV', method= 'oscorespls')


#predict new dN from Landsat 8

#just L8 retrievals - all pixels
L8_ret <- filter(data, SAT == 'L8')
L8_ret_bands <- select(L8_ret, coastal, blue, green, red, NIR, SWIR1, SWIR2)/10000

pred_fit <- predict(L8_wet_fit, ncomp = 7, newdata = data.matrix(L8_ret_bands))

#combine predictions with input data
L8_comb <- mutate(L8_ret, PRED= as.numeric(unlist(pred_fit)))




# Group predicted values by site and average ------------------------------


#function to calculate average predicted values for each site
avg_fun <- function(INPUT)
{
  OUT_comb <- c()
  sites <- unique(INPUT$CODE)
  
  for (i in 1:length(sites))
  {
    #i <- 8
    
    temp_site <- filter(INPUT, CODE == sites[i])
    scenes <- unique(temp_site$SOURCE)
    
    OUT_comb_scene <- c()
    for (j in 1:length(scenes))
    {
      #j <- 8
      temp_scene <- filter(temp_site, SOURCE == scenes[j])
      temp_scene_pred <- as.numeric(unlist(select(temp_scene, PRED)))
      
      pred_mean <- mean(temp_scene_pred)
      pred_sd <- sd(temp_scene_pred)
      
      temp_select <- select(temp_scene, SAT, CODE, YEAR, DAY,
                            SOURCE, latitude, longitude)[1,]
      
      temp_scene_out <- cbind(temp_select, pred_mean, pred_sd)
      OUT_comb_scene <- rbind(OUT_comb_scene, temp_scene_out)
    }
    
    OUT_comb <- rbind(OUT_comb, OUT_comb_scene)
  }
  
  return(OUT_comb)
}



#Run function

#NA vals for Pred_sd indicate only one pixel was retreived at that site/date
L45_out <- avg_fun(L45_comb)
L7_out <- avg_fun(L7_comb)
L8_out <- avg_fun(L8_comb)


#bind all sensors together
master_bound_pre <- rbind(L45_out, L7_out, L8_out)

