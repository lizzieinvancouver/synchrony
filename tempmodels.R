## Started May-ish 2016 ##
## By Heather ##
## Some edits by Lizzie in late May while working on PP checks ##


rm(list=ls())
#row.names=FALSE
# setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
setwd("~/Documents/git/projects/trophsynch/synchrony")

# libraries
library(rstan)

## HAS TEMP CHANGED?
clim<-read.csv("input/climate3.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
clim$name_env<-with(clim, paste(species, phenophase, sep="_"))
clim<-subset(clim, phenophase!="start" & extra!="10" & envunits=="C" & studyid!="HMK028" & studyid!="HMK029" & studyid!="HMK037") #029,037 because nutrients>phenology; 028- because temperature sum
clim$envvalue<-as.numeric(clim$envvalue)
clim$envfactor[clim$envfactor=="temperaure"] <- "temperature"

# Mean temp across sites for those studies with multiple sites
sites<-subset(clim, studyid=="HMK018" | studyid=="HMK019" | studyid=="HMK023" & site!="tomakomai")
nosites<-subset(clim, site=="tomakomai" | studyid!="HMK018" & studyid!="HMK019" & studyid!="HMK023")
nosites<-nosites[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")]

new<-with(sites, aggregate(envvalue, by=list(studyid, year, species), FUN=mean, na.rm=T)) # across sites
names(new)[1]<-"studyid"; names(new)[2]<-"year"; names(new)[3]<-"species"; names(new)[4]<-"envvalue"
new<-new[order(new$studyid),]
#'all' species: "HMK018" "HMK023" "HMK028" "HMK029" "HMK036" "HMK042" "HMK043" >> now fixed in climate3.csv

sites2<-merge(sites[,c("studyid","envfactor","envunits","envtype","year","species")], new, by=c("studyid","year","species"))
sites3<-unique(sites2[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")])

clim2<-rbind(nosites, sites3) # all years with data
#merge with spp data so calculating env change only over the years of interaction
total3<-read.csv("input/raw_april.csv", header=TRUE, na.strings="NA", as.is=TRUE)
total3<-na.omit(total3)
clim3<-merge(clim2, total3[,c("studyid","year")], by=c("studyid","year"))
clim3<-unique(clim3[,c("studyid","year","envfactor","envunits","envtype","species","envvalue")])
clim3<-na.omit(clim3)

#Hinge
hinge<-subset(clim3, species=="Acrocephalus arundinaceus" | species=="Acrocephalus scirpaceus" | species=="Copepod1 spp." | species=="Copepod2 spp." | species=="Cyclops vicinus"  | species=="Daphnia3 spp." | species=="Diatom3 spp." | species=="Perca fluviatillis" | species=="Phytoplankton1 spp." | species=="Pleurobrachia pileus" | species=="Pleurobrachia_a pileus")

hinge_non<-subset(clim3, species!="Acrocephalus arundinaceus" & species!="Acrocephalus scirpaceus"  & species!="Copepod1 spp." & species!="Copepod2 spp." & species!="Cyclops vicinus"  & species!="Daphnia3 spp." & species!="Diatom3 spp." & species!="Perca fluviatillis" & species!="Phytoplankton1 spp." & species!="Pleurobrachia pileus" & species!="Pleurobrachia_a pileus")

hinge_non$newyear<-hinge_non$year
hinge_pre<-subset(hinge, year<=1981); hinge_pre$newyear<-1981
hinge_post<-subset(hinge, year>1981); hinge_post$newyear<-hinge_post$year
hinges<-rbind(hinge_pre, hinge_post)

clim4<-rbind(hinge_non, hinges);
clim3<-clim4

#Model
clim3$yr1981 <- clim3$newyear-1981
N <- nrow(clim3)
y <- clim3$envvalue
specieschar.hin<- aggregate(clim3["envvalue"], clim3[c("studyid", "species")], FUN=length) #number of years per species
J <- nrow(specieschar.hin)
species <- as.numeric(as.factor(clim3$species))
year <- clim3$yr1981
nVars <-1
Imat <- diag(1, nVars)

#New model as of April 2016
#Random slopes only, no random intercepts, hinge, no covariate matrix
# Note to HMK: I don't believe you need the Imat if you have no covariate matrix, right?
temp.model<-stan("stan/synchrony_apr_nocov.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)


###
### Posterior predictive checks
###

# get the lm fits 
temptrends <- read.csv("output/lmfits_temp.csv", header=TRUE)

# First, plot the real raw data used in the model 
pdf("graphs/ppchecks/tempbytime/realdata_formodel.pdf", height=4, width=6)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw real temp data")
for (j in 1:J){
  lines(year[species==j], y[species==j])
}
hist(y, xlab="Day of year", main="Real temp data, cleaned (I assume)")
dev.off()

# Now grab the stan output 
goo <- extract(temp.model) 
hist(goo$mu_b) # example!

# extract means for now #
sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)
sigma_b <- mean(goo$sigma_b) # from stan output 
mu_b <- mean(goo$mu_b) # from stan output 

a <- colMeans(goo$a) # this may not be ideal but deals with the non-pooling 
b <- rnorm(J, mean=mu_b, sd=sigma_b) # this is one way to create fake data from the Stan output to use in the PP check

# Create the PP data using new a and b for each of 71 species
# This is just creating basically one new set of data, or one predictive check
ypred <- length(N) # Lizzie added
for (n in 1:N){
    s <- species[n]
    ypred[n] <- a[s] + b[s]*year[n] # one way to handle the non-pooled intercepts, there may be other ways
}
y <- rnorm(N, ypred, sigma_y)

# Plot the data and see what the raw data predicted from the model looks like
pdf("graphs/ppchecks/tempbytime/onepredictivecheck.pdf", height=4, width=6)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year",
    bty="l", main="Data from posterior means")
for (j in 1:J)
  lines(year[species==j], y[species==j])
hist(y, xlab="Day of year", main="Data from posterior means")
dev.off()
# Hmm, the model makes more data than we seem to have ...
# but we assume that's because of repeating data issues (some lines in raw data figure are on top of each other)

# compare a few things on this single new dataset
par(mfrow=c(2,2))
hist(a, main="interecepts (a) from the stan model with mean from the raw data in blue")
abline(v = mean(temptrends$intfits), col = "blue", lwd = 2) # not even on the graph! ** Are the YEARS centered differently between the model and the raw data?
hist(temptrends$intfits, main="intercepts from the lm fits")
hist(b, main="slopes (b) from the stan model with mean from the raw data in blue")
abline(v = mean(temptrends$slopefits), col = "blue", lwd = 2) # less negative, slopes are generally pooled towards center which makes sense
hist(temptrends$varfits, main="sigma y (b) from the raw data with sigma_y from the model in blue")
abline(v=sigma_y, col = "blue", lwd = 2) 

# okay, but predictive checks are about much more than ONE simulation, one draw ...
# Create the data using new a and b for each of 71 species, 100 times
y.sd100 <- matrix(0, ncol=100, nrow=J)
for (i in 1:100){
    for (n in 1:N){
        s <- species[n]
        ypred[n] <- goo$a[s] + b[s]*year[n] # I think a[s] would also work
    }
  y <- rnorm(N, ypred, sigma_y)
  y.df <- as.data.frame(cbind(y, species))
  y.sd <- aggregate(y.df["y"], y.df["species"], FUN=sd)
  y.sd100[,i] <- y.sd[,2] 
}
hist(colMeans(y.sd100))
# and here's the real data
real.sd <- aggregate(clim["envvalue"], clim[c("studyid", "species")],
    FUN=sd)

hist(colMeans(y.sd100), col=rgb(1,0,0,0.5),
    main="Mean(SD of y) from 100 simulated datasets based on Stan model",
    ylab="mean(SD of one sim)")
abline(v = mean(real.sd$envvalue, na.rm=TRUE), col = "blue", lwd = 2)

## looking at slopes
b100 <- matrix(0, ncol=100, nrow=J)
for (j in 1:100){
    b100[,j] <- rnorm(J, mean=mu_b, sd=sigma_b)
}

hist(colMeans(b100), col=rgb(1,0,0,0.5), xlim=c(-0.1, 0.2), 
    main="Mean(b) from 100 random draws based on Stan model",
    ylab="mean(b from one draw)")
abline(v = mean(temptrends$slopefits), col = "blue", lwd = 4) # mean of lm fits
hist(temptrends$slopefits, col=rgb(0,0,1,0.5), add=TRUE)
# Some very strong pooling of the slopes it looks like, both the negative fits and the positive high fits
# Is this what we want? Might be good to talk to Stan group about this
# Also we should check that the pooling appears correct (e.g., are the shortest time series being pooled the most?)

###
### end posterior predictive checks
###

##
#Exploratory figures
ggplot(subset(clim3, envtype=="air"), aes(x=year, y=envvalue))+
geom_point(aes(colour=factor(species)))+xlim(1969, 2013)+geom_vline(xintercept=1981, colour="grey", linetype = "longdash")+theme_bw()#+geom_smooth(method="lm", se=FALSE, aes(colour = factor(species)))

#by species
ggplot(subset(clim3, envtype=="air"), aes(x=year, y=envvalue))+
geom_point()+xlim(1969, 2013)+geom_vline(xintercept=1981, colour="grey", linetype = "longdash")+facet_wrap(~species)+theme_bw()

ggplot(subset(clim3, envtype=="water"), aes(x=year, y=envvalue))+
geom_point()+xlim(1969, 2013)+geom_vline(xintercept=1981, colour="grey", linetype = "longdash")+facet_wrap(~species)+theme_bw()

#by unique weather station
#air- 
ggplot(subset(clim3, envtype=="air"), aes(x=year, y=envvalue))+
geom_point(aes(colour=factor(studyid)))+facet_wrap(~sitecode)
#water- by unique weather station
ggplot(subset(clim3, envtype=="water"), aes(x=year, y=envvalue))+
geom_point(aes(colour=factor(species)))+facet_wrap(~sitecode)
#geom_path(data=yo, aes(x=SiteTypeBinary, y=log(std), group=factor(sppid), colour=factor(sppid)), size=1)+
#geom_point(data=yo, aes(x=SiteTypeBinary, y=log(std), colour=factor(sppid)))

potential hinge
air:
Acrocephalus arundinaceus (207)
Acrocephalus scirpaceus (207)
water:
Copepod1 spp. (193)
Copepod2 spp. (208)
Cyclops vicinus (191)
Daphnia3 spp.
Diatom3 spp.
Perca fluviatillis
Phytoplankton1 spp.
Pleurobrachia pileus
Pleurobrachia_a pileus
