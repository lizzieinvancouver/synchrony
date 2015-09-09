## Started 20 August 2015 ##
## By Lizzie ##

## Trying to do predictive checks from Stan output ##

setwd("~/Documents/git/projects/trophsynch/synchrony")
source("syncmodels.R")
setwd("~/Documents/git/projects/trophsynch/synchrony")
source("syncmodels.R")

# First, plot the real data used in the model (no dups)
pdf("graphs/realdata_formodel.pdf", height=4, width=6)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw real data (after cleaning)")
for (j in 1:J){
  lines(year[species==j], y[species==j])
}
hist(y, xlab="Day of year", main="Real data, cleaned")
dev.off()

# look at output 
goo <- extract(fit.notypewcov)
hist(goo$mu_b) # example!

# extract means for now #
sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)
sigma_a <- mean(goo$sigma_a) # from stan output 
sigma_b <- mean(goo$sigma_b) # from stan output 
mu_b <- mean(goo$mu_b) # from stan output 
mu_a <- mean(goo$mu_a) # from stan output 

a <- rnorm(J, mean=mu_a, sd=sigma_a) # alert! perhaps should not set sd to sigma exactly?
b <- rnorm(J, mean=mu_b, sd=sigma_b) # alert! perhaps should not set sd to sigma exactly?

# Create the data using new a and b for each of 71 species
ypred <- length(N) # Lizzie added
for (n in 1:N){
    s <- species[n]
    ypred[n] <- a[s] + b[s]*year[n]
}
y <- rnorm(N, ypred, sigma_y)

# Plot the data
pdf("graphs/onepredictivecheck.pdf", height=4, width=6)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year",
    bty="l", main="Data from posterior means")
for (j in 1:J)
  lines(year[species==j], y[species==j])
hist(y, xlab="Day of year", main="Data from posterior means")
dev.off()


## The above uses the year, N and J of the *real data* ##
## The below would let you change the year but I found it unrealistically gave
# too many longer time series so for now just using N, J and year from the real data ##
year_0 <- 1976 
n_data_per_species <- round(runif(J, 6, 37)) 
species <- rep(1:J, n_data_per_species)
N <- length(species)
year <- rep(NA, N)
for (j in 1:J){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0
}

