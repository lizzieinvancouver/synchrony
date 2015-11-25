## Started 20 August 2015 ##
## By Lizzie ##

## Trying to do predictive checks from Stan output ##

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

# get the lm fits
spptrends <- read.csv("output/synchrony1spptrends.csv", header=TRUE)

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

## Last we met I was told to ....
# mean of species' slopes from real data versus slopes from 100 runs of simulated data
# sd of y versus sd of y from the 100 runs of simulated data
## Add in uniform distribution for a for stan model and see if that improves anything (tried, won't run)


# Create the data using new a and b for each of 71 species, 100 times
y.sd100 <- matrix(0, ncol=100, nrow=J)
for (i in 1:100){
    for (n in 1:N){
        s <- species[n]
        ypred[n] <- a[s] + b[s]*year[n]
    }
  y <- rnorm(N, ypred, sigma_y)
  y.df <- as.data.frame(cbind(y, species))
  y.sd <- aggregate(y.df["y"], y.df["species"], FUN=sd)
  y.sd100[,i] <- y.sd[,2] 
}
hist(colMeans(y.sd100))
# and here's the real data
real.sd <- aggregate(rawlong.nodups["phenovalue"], rawlong.nodups[c("studyid", "spp")],
    FUN=sd)

hist(colMeans(y.sd100), col=rgb(1,0,0,0.5), xlim=c(11.2,13), ylim=c(0,35),
    main="Mean(SD of y) from 100 simulated datasets based on Stan model",
    ylab="mean(SD of one sim)")
abline(v = mean(real.sd$phenovalue), col = "blue", lwd = 2)
# Overlay option ... need to adjust X to use
# hist(real.sd$phenovalue, col=rgb(0,0,1,0.5), add=TRUE)

b100 <- matrix(0, ncol=100, nrow=J)
for (j in 1:100){
    b100[,j] <- rnorm(J, mean=mu_b, sd=sigma_b)
}
# get the lm fits
spptrends <- read.csv("output/synchrony1spptrends.csv", header=TRUE)

hist(colMeans(b100), col=rgb(1,0,0,0.5), xlim=c(-0.7, -0.1), ylim=c(0,35),
    main="Mean(b) from 100 random draws based on Stan model",
    ylab="mean(b from one draw)")
abline(v = mean(spptrends$lmfits), col = "blue", lwd = 4)
# Overlay option ... need to adjust xlim to c(-3, 1.5) 
# hist(spptrends$lmfits, col=rgb(0,0,1,0.5), add=TRUE)


hist(apply(b100, 2, sd), col=rgb(1,0,0,0.5), xlim=c(0.4, 0.7), ylim=c(0,30),
    main="SD(b) from 100 random draws based on Stan model",
    ylab="SD(b from one draw)")
abline(v = sd(spptrends$lmfits), col = "blue", lwd = 4)


####################
####################

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

