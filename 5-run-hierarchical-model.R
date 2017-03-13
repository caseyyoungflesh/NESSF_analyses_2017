######################
#Characterizing change in prediced dN over time at Adelie penguin colonies
#
#Script to run model for change in predicted dN over time
#
#Authors: Casey Youngflesh
######################



# Process data for JAGS ---------------------------------------------------

#assign id # for codes to run through in JAGS model
codes <- unique(master_bound_pre$CODE)
cd_df <- data.frame(CODE = codes, id= 1:length(codes))
master_bound_pre2 <- merge(master_bound_pre, cd_df, by = 'CODE')
master_pred <- select(master_bound_pre2, SAT, CODE, S_ID = id, 
                      YEAR, DAY, SOURCE, latitude, longitude, pred_mean, pred_sd)



J_data <- list(
  y = master_pred$pred_mean, #reponse variables
  id = master_pred$S_ID, #column of group id
  N = length(master_pred$pred_mean), #length of data
  M = length(unique(master_pred$CODE)), #number of colonies
  x =  c(1:length(master_pred$pred_mean))) #covariate



# JAGS model to look at changes in pred dN over time at each site ------------------------


#linear model for each site - modeled hierarchically


setwd(paste0(dir, 'Spec_change_over_time_Schwaller_algo/JAGS_files'))

sink("Pred_dN_change.jags")

cat("
    model {
    
    #id is site CODE
    #x is year
    #y is mean predicted dN at a site
    #M is number of ids (sites)
    #N is number of datapoints
    
    for (i in 1:N)
    {
    y[i] ~ dnorm(mu.y[i], tau.y)
    
    #either separate betas for each site (hierarchical), or single beta (and separate alphas)
    mu.y[i] <- alpha[id[i]] + beta[id[i]]*x[i] 
    }
    
    #priors
    
    #for each id
    for (j in 1:M)
    {
    alpha[j] ~ dnorm(mu.a, tau.a)
    beta[j] ~ dnorm(mu.b, tau.b) 
    }
    
    #priors
    tau.y ~ dgamma(.01,.01) #shape and rate
    mu.a ~ dnorm(0, 0.001)
    mu.b ~ dnorm(0, 0.001)
    tau.a ~ dgamma(0.01, 0.01)
    tau.b ~ dgamma(0.01, 0.01)
    
    }",fill = TRUE)
  
sink()



#----------------------#
#Starting values

InitStage <- function() {list(alpha = rep(0.1, length(unique(master_pred$CODE))),
                              beta = rep(0.1, length(unique(master_pred$CODE))),
                              tau.y = 0.1,
                              
                              mu.a = dnorm(1),
                              mu.b = dnorm(1),
                              tau.a = 0.1,
                              tau.b = 0.1
                              
)}



#----------------------#
#Parameters to track

#for single species

ParsStage <- c('alpha', 
               'beta',
               'mu.a',
               'mu.b'
)



#----------------------#
#Inputs for MCMC

ni <- 200000  # number of draws from the posterior
nt <- 2    #thinning rate
nb <- 190000  # number to discard for burn-in 
nc <- 3  # number of chains



#----------------------#
#Run model

out = jags(inits=InitStage,
           n.chains=nc,
           model.file="Pred_dN_change.jags",
           working.directory=getwd(),
           data=J_data, 
           parameters.to.save=ParsStage,
           n.thin=nt,
           n.iter=ni,
           n.burnin=nb,
           DIC=T)



# Analyze model output ----------------------------------------------------

MCMCsummary(out, params = 'beta')
MCMCtrace(out, params = 'beta', pdf = TRUE)
MCMCplot(out, params = 'beta', rank = TRUE, labels = NULL,
         thick_sz = 2, thin_sz = 1, med_sz = 1,
         xlim = c(-0.1, 0.1))

#extract chains for betas - merge the [X] with id
beta_chains <- MCMCchains(out, params = 'beta')
beta_medians <- apply(beta_chains, 2, median)
col_meds <- data.frame(S_ID = 1:max(master_pred$S_ID), beta_medians)

MCMCsummary(out, params = 'mu')
MCMCtrace(out, params = 'mu', pdf = FALSE)
MCMCplot(out, params = 'mu.b', rank = TRUE)




# Prep data for mapping ---------------------------------------------------

OUT <- c()
p_sites <- unique(master_pred$S_ID)
for(i in 1:length(p_sites))
{
  #i <- 1
  temp <- filter(master_pred, S_ID == p_sites[i])
  t_mn <- mean(temp$pred_mean)
  t_cm <- filter(col_meds, S_ID == p_sites[i])[,2]
  
  t_master <- filter(master_pred, S_ID == p_sites[i])[1,]
  t_pm <- select(t_master, CODE, S_ID, latitude, longitude)
  
  t_out <- data.frame(t_pm, mn_dN = t_mn, med_b = t_cm)
  
  OUT <- rbind(OUT, t_out)
}

t_plot <- OUT




# Define map -------------------------------------------------------------


#credit to C. Che-Castaldo for mapping code
world <- map_data('world')
worldmap <- ggplot(world, aes(x = long, y = lat, group = group)) + 
  geom_polygon(color = "black", fill = "white") + 
  geom_path() + 
  scale_y_continuous(breaks = (-2:2) * 30) + 
  scale_x_continuous(breaks = (-4:4) * 45) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank())

antarctica <- worldmap + coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90, -60))



# Plot mean predicted dN for each site ------------------------------------

#compare spatial differences in average predicted dN
antarctica + 
  geom_point(data = t_plot, 
             aes(x = longitude, 
                 y = latitude, 
                 group = S_ID, 
                 colour = mn_dN), 
             size = 2) +
  scale_colour_gradient2("Mean dN", 
                        limits = c(min(t_plot$mn_dN), 
                                   max(t_plot$mn_dN)), 
                        low = "red", 
                        high = "blue") + 
  theme(legend.position = c(0.15, 0.15))




# Plot estimated change for each site -------------------------


antarctica + 
  geom_point(data = t_plot, 
             aes(x = longitude, 
                 y = latitude, 
                 group = S_ID, 
                 colour = med_b), 
             size = 4, alpha = 0.8) +
  scale_colour_gradient2("Estimated Tr", 
                        limits = c(min(t_plot$med_b), 
                                   max(t_plot$med_b)), 
                        low = "red", 
                        high = "blue") + 
  theme(legend.position = c(0.1, 0.15))




# Compare diet trends of AP vs RS -----------------------------------------


AP <- filter(t_plot, longitude > -110 & longitude < -50)
EA <- filter(t_plot, longitude > 0 & longitude < 179)


#differences in dN trend between regions
t.test(AP$med_b, EA$med_b)

mean(AP$med_b) #decreasing on AP
mean(EA$med_b) #increasing in EA


#differences in mean dN between regions
t.test(AP$mn_dN, EA$mn_dN)

mean(AP$mn_dN) #lower in WAP
mean(EA$mn_dN) #higher in East Antarctica


