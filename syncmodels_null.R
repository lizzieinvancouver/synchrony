## Started 21 August 2015 ##
## By Lizzie ##

## Trying to do a quick null model from Stan output ##
## Based off current code from syncmodels_PPchecks.R ##

setwd("~/Documents/git/projects/trophsynch/synchrony")
source("syncmodels.R")

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
m <- matrix(0, ncol = 20, nrow = nrow(rawlong.nodups))
ypreds <- data.frame(m)
ypred <- length(N)
for (x in c(1:ncol(ypreds))){
    for (n in 1:N){
        s <- species[n]
        ypred[n] <- a[s] + b[s]*year[n]
    }
    y <- rnorm(N, ypred, sigma_y)
    ypreds[,x] <- y
}
quicknull <- cbind(rawlong.nodups, ypreds)

