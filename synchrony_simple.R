# Simulate fake data

# setwd("/Users/Lizzie/Documents/git/projects/trophsynch")

# Create the species-level parameters
J <- 50 # number of species
type <- rep(c(1,2), c(25,25)) # consumer or resource (1 or 2)
level <- c(125, 125)
sigma_level <- c(25, 25)
species_level <- rnorm(J, 0, 1)

## no trend in simple
## so commenting out below
# trend <- c(-.1, -.2)
# sigma_trend <- c(.2, .2)
# species_trend <- rnorm(J, 0, 1)

year_0 <- 1976
sigma_y <- 5

# Create the data
n_data_per_species <- round(runif(J, 5, 40)) # how many years per sp.?
species <- rep(1:J, n_data_per_species)
N <- length(species)
year <- rep(NA, N)
for (j in 1:J){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0
}
ypred <- length(N) # Lizzie added
for (n in 1:N){
  s <- species[n]
  t <- type[s]
  ypred[n] <- level[t] + sigma_level[t]*species_level[s]
}
y <- rnorm(N, ypred, sigma_y);

# Plot the data

pdf("raw1simple.pdf", height=4, width=6)
colors=c("blue", "red")
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw data (types 1 and 2 are blue and red)")
for (j in 1:J)
  lines(year[species==j], y[species==j], col=colors[type[j]])
dev.off()

# Fit the model

library("rstan")
source("stan.R")

fit_simple <- stan("synchrony_simple.stan", data=c("N","y","J","species","year","type"), iter=2000, chains=4)
print(fit_simple)

# Look at the results
library("shinyStan")
launch_shinystan(fit_simple)
