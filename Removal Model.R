
library(tidyverse)
library(rjags)
load.module("mix")
load.module("dic")
library(R2jags)
library(coda)
library(FSA)

setwd("H:/Removal_Model")
seed <- 111
set.seed(seed)

# import dataset that includes fish pass data
fish.data <- read.csv("removal_data.csv")
head(fish.data)

# evaluate the length-frequency distribution
ggplot(fish.data, aes(x = Length)) +
  geom_histogram(binwidth = 20)
ggplot(fish.data, aes(x = log2(Length))) +
  geom_histogram(binwidth = 0.5)

# log transform length
fish.data$log2length <- log2(fish.data$Length)

# assign each fish to a size class
fish.data$lencat[fish.data$log2length < 5 & fish.data$log2length >= 4] <- 4.5
fish.data$lencat[fish.data$log2length < 6 & fish.data$log2length >= 5] <- 5.5
fish.data$lencat[fish.data$log2length < 7 & fish.data$log2length >= 6] <- 6.5
fish.data$lencat[fish.data$log2length < 8 & fish.data$log2length >= 7] <- 7.5
# remove any individuals not assigned a size class
fish.data <- fish.data[!is.na(fish.data$lencat),]

# summarize the fish data by site and size class
data <- fish.data %>%
  dplyr::group_by(Site, lencat) %>%
  summarize(TotalCatch = n(),
            Pass1 = sum(Pass == 1),
            Pass2 = sum(Pass == 2),
            Pass3 = sum(Pass == 3),
            Weight.geomean = geomean(Weight, na.rm = TRUE)) %>%
  droplevels()  %>%
  ungroup() %>%
  complete(Site, lencat) # make sure to include every possible site-size class combination

# number of sites
Nsite <- length(unique(data$Site))
# number of size classes
Nsizeclass <- length(unique(data$lencat))

# Reformat pass catches into matrix form
c1 <- structure(.Data = c(data$Pass1), .Dim = c(Nsite, Nsizeclass))
c2 <- structure(.Data = c(data$Pass2), .Dim = c(Nsite, Nsizeclass))
c3 <- structure(.Data = c(data$Pass3), .Dim = c(Nsite, Nsizeclass))
# geometric mean weight matrix to estimate biomass
gmw <- structure(.Data = c(data$Weight.geomean), .Dim = c(Nsite, Nsizeclass))
gmw <- ifelse(is.na(gmw), 0, gmw) # replace NAs with 0

# values to represent the different size classes
(size <- c(unique(data$lencat)))

# JAGS sampler parameters
nc <- 3 # number of chains
ni <- 200000 # number of iterations/chain length
nb <- 50000 # burn in value
nt <- 20 # thinning rate

# model file
sink("removal.model.bug")
cat("model {
    
    for (i in 1:Nsite) {
      for (j in 1:Nsizeclass) {
    
      c1[i,j] ~ dbin(q[i,j], n[i,j]) # pass 1 catch
      n2[i,j] <- n[i,j] - c1[i,j] # remaining population after pass 1
      c2[i,j] ~ dbin(q[i,j], n2[i,j]) # pass 2 catch
      n3[i,j] <- n2[i,j] - c2[i,j] # remaining population after pass 2
      c3[i,j] ~ dbin(q[i,j], n3[i,j]) # pass 3 catch
  
      n[i,j] <- round(nz[i,j])
      nz[i,j] ~ dlnorm(mu.n[j], tau.n[j]) # initial population size

      # capture probability function
      logit(q[i,j]) <- b0[i] + b.size*size[j] 

      bio[i,j] <- n[i,j] * gmw[i,j] # size class biomass
      
      # predicted values
      c1_pred[i,j] <- n[i,j] * q[i,j]
      c2_pred[i,j] <- n2[i,j] * q[i,j]
      c3_pred[i,j] <- n3[i,j] * q[i,j]
      }

    # random site effect
    b0[i] ~ dnorm(b0.mu, b0.tau)

    # site population size
    nt[i] <- sum(n[i,]) 
    # site biomass
    bio_site[i] <- sum(bio[i,]) 
    
    }
    
    # population size hyperparameters
    for (y in 1:Nsizeclass) {
    mu.n[y] ~ dunif(0,20) # mean
    tau.n[y] ~ dunif(0.0000001, 1) # precision
    sigma.n[y] <- 1/sqrt(tau.n[y]) # standard deviation
    }

    # priors
    b0.mu ~ dnorm(0, 0.0001)
    b0.tau ~ dunif(0,5)
    b.size ~ dnorm(0, 0.0001)

    }")
sink()

# data to input into the model
data <- list("Nsite" = Nsite, "Nsizeclass" = Nsizeclass,
             "c1" = c1, "c2" = c2, "c3" = c3, "gmw" = gmw, "size" = size)

# model parameters to be monitored
params <- c("n", "nt","mu.n", "q", "sigma.n","b.size", "b0", "b0.mu", "b0.tau", "bio", 
            "bio_site", "c1_pred", "c2_pred", "c3_pred") 

# save model output as an object
jagsfit <- R2jags::jags(data = data, parameters.to.save = params, n.chains = nc,
                        n.iter = ni, n.thin = nt, n.burnin = nb, 
                        model.file = "removal.model.bug")

print(jagsfit)
