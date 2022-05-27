library(jagsUI)

rm(list=ls())

setwd("~/iles_ECCC/Landbirds/Boreal-monitoring/YT-GHCH-COSEWIC/GHCH_population_size")

# ------------------------------------------------
# Analysis method 2: 
#  - Use informative prior, based on estimates of density from elsewhere
# ------------------------------------------------

# Model script
sink("ghch_model.jags")
cat("

    model {

      # For informative prior
      log_median <- log(0.1)
      log_sigma <- 1.5
      log_tau <- pow(log_sigma,-2)
      male_dens ~ dlnorm(log_median,log_tau)
      
      total_area_detectable <- IAO_area * male_dens * pi * (territory_radius + EDR)^2
      total_area_occupied <- IAO_area * male_dens * pi * (territory_radius)^2
      proportion_detectable <- total_area_detectable/IAO_area
      
      # Likelihood
      y ~ dbin(proportion_detectable, nARU)
      
      # Total population size of mature individuals
      popsize <- IAO_area * male_dens * 2
    }
    
",fill = TRUE)
sink()

bird_EDR = 0.04            # km
bird_Territory_Area = 0.17 # km2
bird_Territory_Radius = sqrt(bird_Territory_Area/pi)
bird_Territory_Radius_detectable = bird_Territory_Radius + bird_EDR
nARU = 82

area_monitored = AD*nARU # Total area monitored by all ARUs
IAO_area = 3303

# Analyze using JAGS
jags.data <- list(y = 0,   # Number of ARUs that detected the species
                  nARU = nARU,
                  IAO_area = IAO_area,
                  territory_radius = bird_Territory_Radius,
                  EDR = bird_EDR,
                  pi = pi)

out <- jags(data = jags.data,
            model.file = "ghch_model.jags",
            parameters.to.save = c("male_dens","popsize"),n.chains = 3,
            n.thin = 1,
            n.iter = 210000,
            n.burnin = 10000)

# Examine results
hist(out$sims.list$male_dens, breaks = seq(0,max(out$sims.list$male_dens)+1,0.005), xlim = c(0,1), border = "transparent",
     main = "Posterior estimate of density (males / km2)", xlab = "Density (males / km2)")
mean(out$sims.list$male_dens)

# This is the estimate of total number of mature individuals in IAO
hist(out$sims.list$popsize, breaks = seq(0,max(out$sims.list$popsize)+10,10), xlim = c(0,2000), border = "transparent",
     main  = "Posterior estimate of abundance (males)", xlab = "Abundance of males")
median(out$sims.list$popsize)

# Calculate of total number of individuals
p_less_than_1000 <- mean(out$sims.list$popsize < 1000)
p_less_than_1000

p_less_than_250 <- mean(out$sims.list$popsize < 250)
p_less_than_250

quantile(out$sims.list$popsize,c(0.025,0.5,0.975))

