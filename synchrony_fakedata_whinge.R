# Simulate fake data

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")

# RANDOM SLOPE AND INTERCEPT MODEL #
# Create the species-level parameters
J <- 30 # number of species 
#level <- 120 # mean phenovalue and SE is 39
a <- rnorm(J, 120, 41) #species error, intercept
#trend <- rep(-0.3,J) # mean slope?; range is -1.75 to 1.21 chosen range based on slopes from synchrony_notype_hk
#sigma_trend <- rep(0.005,J)  # variance around slope? range is 0.0019-0.0096 chosen range from SE around slope from synchrony_notype_hk
b <- rnorm(J, -0.32, 0.499) #n=87 species
year_0 <- 1981 #nohinge
sigma_y<- 5 #??

# Create the data
n_data_per_species <- round(runif(J, 5, 37)) #how many years per spp
species <- rep(1:J, n_data_per_species) #nrow of dataframe
N <- length(species)
year <- rep(NA, N)
for (j in 1:J){
  year[species==j] <- rev(2013 - 1:(n_data_per_species[j]))-year_0 #from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}


## NO HINGE ##
ypred <- length(N) # Lizzie added
for (n in 1:N){ #actual STAN model
  s <- species[n] #sppid for each row
  ypred[n] <- a[s] + b[s]*year[n]
   }
y <- rnorm(N, ypred, sigma_y);

pdf("raw1_nohinge.pdf", height=4, width=6)
colors=c("blue", "red")
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw data"); abline(b=0)
for (j in 1:J)
  lines(year[species==j], y[species==j])
dev.off()

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
  s <- species[n] #sppid for each row
  ypred[n] <- a[s] + b[s]*year[n]
   }
y <- rnorm(N, ypred, sigma_y);

# Plot the data
pdf("raw1_hinge.pdf", height=4, width=6)
colors=c("blue", "red")
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(x=year, y=y, type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw data", xlim=c(-5,31), ylim=c(39,209));
for (j in 1:J)
  lines(year[species==j], y[species==j])
dev.off()

