### Started 1 April 2015 ###
### April Fool's Day! ###

## Trying to run STAN with synchrony data ##

## Last updated on Tax Dax in the USA ##

options(stringsAsFactors = FALSE)

setwd("~/Documents/git/projects/trophsynch/synchrony")
library(ggplot2)
library(rstan)
# library(lme4)

# get data
# in the data file: spp1 = neg= negative species e.g. resource 
raw <- read.csv("output/raw.csv", header=FALSE)
mat.mm <- read.csv("output/alldata_long.csv", header=TRUE)

##
## lme4: fails to converge
# modresource <- lmer(neg_phenovalue ~ year + (1 + year|spp1), data=raw)
# modresource <- lmer(neg_phenovalue ~ year + (year|spp1), data=raw)

##
## try Stan #
# 3 stages of cooking: prepared foods, recipes, make your own
# Stan does not currently allow ragged arrays

# prep the data to fit the model including:
# aggregate to get species level characteristics
# subset down to the phenovalues
longpheno <- subset(mat.mm, whatmeasured=="phenovalue")
longpheno$yr1976 <- longpheno$year-1976
specieschar <- aggregate(longpheno["measurement"],
    longpheno[c("studyid", "value", "variable")], FUN=length)
specieschar.nodup <- specieschar[!(duplicated(specieschar$value)==TRUE),] # there are some species with multiple relationships

N <- nrow(longpheno)
y <- longpheno$measurement
J <- nrow(specieschar.nodup)
species <- as.numeric(as.factor(longpheno$value))
year <- longpheno$yr1976
type <- as.numeric(as.factor(specieschar.nodup$variable)) 

# Run the model
# source("source/stan.R") # fails for Lizzie: all scheduled cores encountered errors in user code
fit1 <- stan("stan/synchrony1.stan", data=c("N","y","J","species","year","type"), iter=2000, chains=4)

qq<-summary(fit1)
class(qq) # gives the overall summary and then each of the four chains (I think)
qq[[1]][,1]

library("shinyStan")
launch_shinystan(fit1)


