### Started 7 June 2016 ###
### By Lizzie 

### Fake data code for synchrony models ###
### A good dose of this taken from synchrony_fakedata_whinge.R ###

## housekeeping
rm(list=ls())
setwd("~/Documents/git/projects/trophsynch/synchrony/testing123")
library(rstan)
# library(shinystan)

##
## starting simple
## hinge model with *one* hierarchical level

# RANDOM SLOPE AND INTERCEPT MODEL #
# Create the species-level parameters

# Important note about the below -- when the errors are set very high
# such as is true below the model will struggle to get too close to the mu 
# values given; so...
# best to check output with lower sigmas, as well as closer to reality
# also consider running model longer to get the n_eff values up higher.

J <- 50 # number of species 
a <- rnorm(J, 120, 41) #species error, intercept
b <- rnorm(J, -0.35, 0.499) # n=87 species
year_0 <- 1981 # hinge
sigma_y <- 42 #

# Create the data
n_data_per_species <- round(runif(J, 5, 37)) #how many years per spp
species <- rep(1:J, n_data_per_species) #nrow of dataframe
N <- length(species)
year <- rep(NA, N)
for (j in 1:J){
  year[species==j] <- rev(2013 - 1:(n_data_per_species[j]))-year_0 
}

### HINGE MODEL ###
for (j in 1:J){
  w<-year[species==j]<=0 # or 1981
  x<-which(w==FALSE)
  k<-length(w)-length(x)
  if (k>=5){
  	d<-year[species==j]
  	d[1:k]<-0 # or 1981
  	year[species==j]<-d
  }
  }

ypred <- length(N) # Lizzie added
for (n in 1:N){ #actual STAN model
  s <- species[n] # sppid for each row
  ypred[n] <- a[s] + b[s]*year[n]
   }
y <- rnorm(N, ypred, sigma_y);


## check the model works and returns the slopes it was given
# first, add in two things you need to run the model
nVars <-1
Imat <- diag(1, nVars)
# then run it!
fit.hinge <- stan("..//stan/synchrony1_notype_wcovar.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)

print(fit.hinge, pars=c("mu_a","mu_b","sigma_a","sigma_y"))



##
## start of model with two hierarchical levels

H <- 10 # number of classes (or some level above species)
J <- 15 # number of species PER class

Hid <- rep(1:H, each=J) # assign each species to a class
mu_Ha <- rnorm(H, 145, 6) # different intercepts for each class

# now we need to make variation in intercepts among species, but build it off class variation
# this means you're trying to model the variation of species with *regards to class*
# which means (generally smaller) differences since they are the diff from a species' class
# it would be something like the below, but I have run out of time to make the loops
mu_a <- mu_Ha[i] + rnorm(J, 10, 3)

# once you have this (intercepts per species, built around intercepts for class)
# you should be in a better position to build the raw data off the code above
