library(sp)
library(rgeos)
library(sf)
library(tidyverse)
library(jagsUI)

rm(list=ls())

# --------------------------------------------------
# Simulate bird territories
# --------------------------------------------------

# Parameters
bird_dens = 0.05        # birds / km2
bird_EDR = 0.045       # bird EDR
bird_TA = 0.17        # km2 territory size for each bird
bird_Tr = sqrt(bird_TA/pi) # radius of bird territories
bird_TrD = bird_Tr + bird_EDR  # Radius of detectable bird territories

# Study area
plot_area = 3303
expected = plot_area*bird_dens  # Expected number of birds
plot_dim = sqrt(plot_area)
dim <- plot_dim/2

# Bird territories
xy = data.frame(x = runif(expected,-dim,dim),y = runif(expected,-dim,dim)) # Centroid of each bird territory
bird_centroids <- SpatialPoints(coords = xy)

# Territories as polygons
bird_T_poly = gBuffer(bird_centroids, byid = TRUE, width = bird_Tr)
bird_T_sf <- st_as_sf(bird_T_poly)

# Add EDR onto edge of territory
bird_TD_poly = gBuffer(bird_centroids, byid = TRUE, width = bird_TrD)
bird_TD_sf <- st_as_sf(bird_TD_poly)

# Plot
bird_radius_m = round(bird_Tr*1000)
ggplot(bird_T_sf) +
  geom_sf(fill = "transparent", col = "blue")+
  scale_fill_manual(values = c("transparent","orangered"),na.value = "transparent",name="Detection")+
  coord_sf(xlim = c(-dim,dim), ylim = c(-dim,dim))+
  xlab("km")+ylab("km")+
  ggtitle(paste0("Bird Density = ",bird_dens," male / km^2\nBird Territory Radius = ",bird_radius_m," m\nBird EDR = ",bird_EDR*1000," m"))


# --------------------------------------------------
# Simulate ARU placement
# --------------------------------------------------
nARU = 82

# SRS sampling

aru_sf = SpatialPoints(coords = data.frame(x = runif(nARU,-dim,dim), y = runif(nARU,-dim,dim))) %>% st_as_sf()
aru_detections <- sf::st_intersects(bird_TD_sf,aru_sf) %>% as.data.frame()
bird_TD_sf$det = "n"

if (nrow(aru_detections) > 0){
  bird_TD_sf$det[aru_detections[,1]] <- "y"
}

ggplot(bird_T_sf) +
  geom_sf(data = bird_TD_sf,aes(fill = det), col = "red", linetype = 2)+
  geom_sf(fill = "transparent", col = "blue")+
  geom_sf(data = aru_sf, shape = 19)+
  scale_fill_manual(values = c("transparent","orangered"),na.value = "transparent",name="Detection")+
  coord_sf(xlim = c(-dim,dim), ylim = c(-dim,dim))+
  xlab("km")+ylab("km")+
  ggtitle(paste0("Bird Density = ",bird_dens," male / km^2\nBird Territory Radius = ",bird_radius_m," m\nBird EDR = ",bird_EDR*1000," m"))


# --------------------------------------------------
# Calculate number of ARUs with detections
# --------------------------------------------------
nARU_det <- length(unique(aru_detections$col.id))
nARU_det

# ----------------------------------
# Repeated simulations
# ----------------------------------

# Model script
sink("ghch_uninformative.jags")
cat("

    model {

      # For informative prior
      log_median <- log(0.1)
      log_sigma <- 1.5
      log_tau <- pow(log_sigma,-2)
      male_dens ~ dlnorm(log_median,log_tau)
      
      # male_dens ~ dgamma(0.001,0.001)
      total_area_detectable <- IAO_area * male_dens * pi * (territory_radius + EDR)^2
      total_area_occupied <- IAO_area * male_dens * pi * (territory_radius)^2
      proportion_detectable <- total_area_detectable/IAO_area
      
      # Likelihood
      y ~ dbin(proportion_detectable, nARU)
      
    }
    
",fill = TRUE)
sink()

nreps <- 1000
detection_sim_vec <- rep(NA,nreps)
dens_estimate_mean_vec <- dens_estimate_median_vec <- rep(NA,nreps)
dens_se_vec <- rep(NA,nreps)
coverage_vec <- rep(NA,nreps)

for (i in 1:nreps){
  
  nARU_det = 0
  
  # Bird territories
  xy = data.frame(x = runif(expected,-dim,dim),y = runif(expected,-dim,dim)) # Centroid of each bird territory
  bird_centroids <- SpatialPoints(coords = xy)
  
  # Territories as polygons
  bird_T_poly = gBuffer(bird_centroids, byid = TRUE, width = bird_Tr)
  bird_T_sf <- st_as_sf(bird_T_poly)
  
  # Add EDR onto edge of territory
  bird_TD_poly = gBuffer(bird_centroids, byid = TRUE, width = bird_TrD)
  bird_TD_sf <- st_as_sf(bird_TD_poly)
  
  aru_sf = SpatialPoints(coords = data.frame(x = runif(nARU,-dim,dim), y = runif(nARU,-dim,dim))) %>% st_as_sf()
  aru_detections <- sf::st_intersects(bird_TD_sf,aru_sf) %>% as.data.frame()

  nARU_det <- length(unique(aru_detections$col.id))
  detection_sim_vec[i] <- nARU_det
  
  # Analyze using JAGS
  jags.data <- list(y = nARU_det,
                    nARU = nARU,
                    IAO_area = plot_area,
                    territory_radius = bird_Tr,
                    EDR = bird_EDR,
                    pi = pi)
  
  out <- jags(data = jags.data,
              model.file = "ghch_uninformative.jags",
              parameters.to.save = c("male_dens"),n.chains = 3,
              n.thin = 1,
              n.iter = 10000,
              n.burnin = 5000)
  
  dens_estimate_mean_vec[i] <- out$mean$male_dens
  dens_estimate_median_vec[i] <- out$q50$male_dens
  dens_se_vec[i] <- out$sd$male_dens
  
  coverage_vec[i] <- out$q2.5$male_dens < bird_dens &  out$q97.5$male_dens > bird_dens
  
  print(i)
}

# Expected number of birds detected
total_area_occupied <- plot_area * bird_dens * pi*bird_TrD^2
proportion_occupied <- total_area_occupied/plot_area
expected_ARU_detections <- nARU*proportion_occupied

# Average number of birds detected across 1000 simulations
mean(detection_sim_vec, na.rm = TRUE) 
expected_ARU_detections # Compare to above

# Model-based estimate (median)
mean(dens_estimate_median_vec, na.rm = TRUE)

# Distribution of bias in estimates
bias_mean <- dens_estimate_mean_vec - bird_dens
bias_median <- dens_estimate_median_vec - bird_dens

mean(bias_mean)
mean(bias_median)

# Model-based 95% interval coverage 
mean(coverage_vec,na.rm = TRUE)

# Model-based SE
mean(dens_se_vec,na.rm = TRUE)
sd(dens_estimate_vec,na.rm = TRUE) # Actual SE

# ----------------------------------
# Analyze using new analysis method
# ----------------------------------

nreps <- 10000
sim_results <- data.frame()

for (i in 1:nreps){
  nsites <- 200
  n_det_site <- rep(NA,nsites)
  
  for (s in 1:nsites){
    n_detect = 0
    
    # Simulate data for this run
    xy = data.frame(x = runif(expected,-dim,dim),y = runif(expected,-dim,dim)) # Centroid of each bird territory
    if(nrow(xy)>0){
      bird_centroids <- SpatialPoints(coords = xy)
      
      bird_TD_poly = gBuffer(bird_centroids, byid = TRUE, width = bird_TrD)
      bird_TD_sf <- st_as_sf(bird_TD_poly)
      
      # ARU placement
      aru_sf = SpatialPoints(coords = data.frame(x = 0, y = 0)) %>% st_as_sf()
      aru_detections <- sf::st_intersects(bird_TD_sf,aru_sf) %>% as.data.frame()
      
      n_detect <- nrow(aru_detections)
    }
    
    n_det_site[s] <- n_detect
  }
  
  jags.data <- list(y = sum(n_det_site),
                    nsites = nsites,
                    AD = pi * (bird_Tr + bird_EDR)^2)
  
  out <- jags(data = jags.data,
              model.file = "new.jags",
              parameters.to.save = c("dens"),n.chains = 3,n.thin = 1,n.iter = 10000,n.burnin = 5000)
  
  df = data.frame(mean = mean(out$sims.list$dens),
                  med = median(out$sims.list$dens),
                  lcl = quantile(out$sims.list$dens,0.025),
                  ucl = quantile(out$sims.list$dens,0.975))
  sim_results <- rbind(sim_results, df)
  
  plot(sim_results$mean ~ seq(1,nrow(sim_results)), type = "p",ylab = "mean density estimate",xlab = "nsims", col = "red", pch = 19)
  abline(h = mean(sim_results$mean), lty = 2, col = "red")
  abline(h = bird_dens,lty = 1, col = "blue")
}

# Average estimate of density
mean(sim_results$mean)

# 95% interval coverage
mean(sim_results$lcl < bird_dens & sim_results$ucl > bird_dens )
