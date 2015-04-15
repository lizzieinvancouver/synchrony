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
raw <- read.csv("output/raw.csv", header=TRUE)
rawlong <- read.csv("output/alldata_long.csv", header=TRUE)

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
longpheno <- subset(rawlong, whatmeasured=="phenovalue")
longpheno$yr1976 <- longpheno$year-1976
specieschar <- aggregate(longpheno["measurement"],
    longpheno[c("studyid", "value", "variable", "intid")], FUN=length)
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

library("shinyStan")
launch_shinystan(fit1)

# Hack to look at the data
qq<-summary(fit1)
class(qq) # gives the overall summary and then each of the four chains (I think)
qq[[1]][,1]
spptrends.fromstan <- as.vector(qq[[1]][,1])[4:73]
spptrends <- cbind(specieschar.nodup, spptrends.fromstan)
spptrends$doybydec <- spptrends$spptrends.fromstan*10
  
ggplot(data=spptrends, aes(x=variable, y=doybydec))+
geom_path(aes(group=factor(intid)))+
    #geom_point(aes(colour=factor(interaction), size=1))+ 
    theme_bw()+theme(panel.grid.major = element_blank(),
    axis.title.x =element_text(size=12), axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, angle=90))+
    ylab("Phenological shift (days/decade)")+xlab("Resource (1) or consumer (2)")+
    ggtitle("Slopes from Stan \n resources are on left, consumers on right (I think)")
