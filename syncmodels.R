### Started 1 April 2015 ###
### April Fool's Day! ###

## Trying to run STAN with synchrony data ##

## Last updated 20 August 2015 ##


# Note that Heather calls synchrony1_notype_wcovar.stan 
# synchrony1_notype_hk.stan

rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/git/projects/trophsynch/synchrony")
library(ggplot2)
library(rstan)
# library(lme4)

# get data
# in the data file: spp1 = neg= negative species e.g. resource 
raw <- read.csv("output/raw.csv", header=TRUE)
rawlong <- read.csv("input/indiv_sppdata2.csv", header=TRUE)
rawlong$int_type[which(rawlong$int_type=="poll")] <- "pollination"
rawlong.X <- NULL

##
# figure out how many unique species
# and alter the data so each unique species shows up once
# watch out: if you leave in int_type you still get duplicated species...
# so it's not below until after duplicate species are removed!
specieschar.wdups <- aggregate(rawlong["phenovalue"],
    rawlong[c("studyid", "species", "spp")], FUN=length)
specieschar <- aggregate(specieschar.wdups["phenovalue"],
    specieschar.wdups[c("studyid", "species")], FUN=length)
dupspp <- subset(specieschar, phenovalue>1)
specieschar.wdups[which(specieschar.wdups$species %in% dupspp$species),]
# delete duplicate species as spp1 (generally the shorter timeseries)
rawlong.nodups <- rawlong[-(which(rawlong$species %in% dupspp$species &
    rawlong$spp=="spp1")),]
# and order it!
rawlong.nodups <- rawlong.nodups[with(rawlong.nodups, order(species, year)),]
specieschar.formodel <- aggregate(rawlong.nodups["phenovalue"],
    rawlong.nodups[c("studyid", "species", "int_type", "spp")], FUN=length)


##
## lme4: fails to converge
# modresource <- lmer(neg_phenovalue ~ year + (1 + year|spp1), data=raw)
# modresource <- lmer(neg_phenovalue ~ year + (year|spp1), data=raw)

##
# get estimate from no pooling (fixed effects) and complete pooling
comppool <- lm(phenovalue~year, data=rawlong.nodups)
comppool.wtype <- lm(phenovalue~year*spp-1, data=rawlong.nodups)
nopool <- lm(phenovalue~year*species, data=rawlong.nodups)
uniquespp <- unique(rawlong.nodups$species)
lmfits <- rep(NA, length(uniquespp))
for (eachsp in 1:length(uniquespp)){
    lmhere <- lm(phenovalue~year, data=subset(rawlong.nodups, species==uniquespp[eachsp]))
    lmfits[eachsp] <- coef(lmhere)[2]
  }

##
# looking at hinges
yearzero <- 1981
rawlong.nodups2<-subset(rawlong.nodups, year <= yearzero & spp=="spp2")
rawlong.nodups2$count<-1
wus <- aggregate(rawlong.nodups2["count"], rawlong.nodups2[c("intid")], FUN=sum)
preccstudies <- subset(wus, count>4)

excludestudies <- c("HMK029", "HMK032", "HMK037")
rawlong.excl <- rawlong[which(!rawlong$studyid %in% excludestudies),]

hinges <- rawlong.excl[which(rawlong.excl$intid %in% preccstudies$intid),]
hinges$newyear <- hinges$year
hinges$newyear[hinges$newyear <= yearzero] <- yearzero

# loop!
# we'll do this for each intid, but since that will replicate species I will feed
# the loop the select the species


makehingeplot <- function(intidhere, specieslist){
    for (i in seq_along(specieslist)){
        subby.nohinge <- subset(hinges, intid==intidhere & species==specieslist[i])
        subby.hinge <- subset(hinges, intid==intidhere & species==specieslist[i])
        pdf(paste("graphs/hingetimeseries/", "hingeat", yearzero, " ", intidhere, 
            specieslist[i], ".pdf", sep=""), height=4, width=6)
        par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
        plot(phenovalue~year, data=subby.nohinge, type="l")
        points(phenovalue~year, data=subby.nohinge, cex=0.6)
        abline((lm(phenovalue~year, data=subby.nohinge)))
        lines(phenovalue~newyear, data=subby.hinge, col="red")
        abline((lm(phenovalue~newyear, data=subby.hinge)), col="red")
        dev.off()
    }
}

makehingeplot("170", c("Phytoplankton1 spp.", "Daphnia3 spp."))
makehingeplot("171", c("Perca fluviatillis"))
makehingeplot("177", c("Caterpillar spp.", "Parus1 major"))
makehingeplot("178", c("Ficedula1 albicollis"))
makehingeplot("179", c("Quercus2 robur", "Caterpillar spp."))
makehingeplot("180", c("Operophtera brumata", "Parus2 major"))

# 1981 only -- add in 191, 193, 201, 207, 208
print("Stop and think, the below plots only work for hinges at/after 1981!")
makehingeplot("191", c("Phytoplankton2 spp.", "Cyclops vicinus"))
makehingeplot("193", c("Diatom3 spp.", "Copepod1 spp."))
makehingeplot("201", c("Guillemots spp.", "Rissa1 tridactyla"))
makehingeplot("207", c("Acrocephalus arundinaceus", "Acrocephalus scirpaceus"))
makehingeplot("208", c("Copepod2 spp.", "Pleurobrachia pileus"))

# get the temperature data
tempdat <- read.csv("input/tempsens_nov3.csv", header=TRUE)
tempdat$newyear <- tempdat$year
tempdat$newyear[tempdat$newyear <= yearzero] <- yearzero
tempdat.sm <- subset(tempdat, select=c("studyid", "year", "envvalue", "z", "newyear"))
tempdat.sm <- tempdat.sm[!duplicated(tempdat.sm), ] # remove duplicated rows
tempstudies <- unique(tempdat.sm$studyid)

for (i in seq_along(tempstudies)){
    subby <- subset(tempdat.sm, studyid==tempstudies[i])
    pdf(paste("graphs/hingetimeseries/", "hingeat", yearzero, " ", tempstudies[i],
        ".pdf", sep=""), height=4, width=6)
    par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
    plot(envvalue~year, data=subby, type="l")
    points(envvalue~year, data=subby, cex=0.6)
    abline((lm(envvalue~year, data=subby)))
    lines(envvalue~newyear, data=subby, col="red")
    abline((lm(envvalue~newyear, data=subby)), col="red")
    dev.off()
}


# end hinges
##


##
## try Stan #hinges170daph3 <- subset(hinges, intid=="170" & species=="Daphnia3")

# 3 stages of cooking: prepared foods, recipes, make your own
# Stan does not currently allow ragged arrays


# prep the data to fit the model including:
# aggregate to get species level characteristics
# subset down to the phenovalues
rawlong.nodups$yr1976 <- rawlong.nodups$year-1976

N <- nrow(rawlong.nodups)
y <- rawlong.nodups$phenovalue
J <- nrow(specieschar.formodel)
species <- as.numeric(as.factor(rawlong.nodups$species))
year <- rawlong.nodups$yr1976
type <- as.numeric(as.factor(specieschar.formodel$spp)) 

# Run the model
# source("source/stan.R")

# Heather's update to Margaret Kosmala's model to add in covar matrix
nVars <- 1
Imat <- diag(1, nVars)
fit.notypewcov <- stan("stan/synchrony1_notype_wcovar.stan", data=c("N","y","J","species","year", "nVars", "Imat"), iter=2000, chains=4)
print(fit.notypewcov)

yy <- summary(fit.notypewcov)
intercepts.fitnotypewcovar.fromstan <- as.vector(yy[[1]][,1])[1:71]
spptrends.fitnotypewcovar.fromstan <- as.vector(yy[[1]][,1])[72:142]

stop("Lizzie made the code stop here")

# Margaret Kosmala's model with no type
fit.notype <- stan("stan/synchrony1_notype.stan", data=c("N","y","J","species","year"), iter=2000, chains=3)
print(fit.notype)

zz<-summary(fit.notype)
intercepts.fitnotype.fromstan <- as.vector(zz[[1]][,1])[1:71]
spptrends.fitnotype.fromstan <- as.vector(zz[[1]][,1])[72:142]

# model with type (resource or consumer); the one and only original
fit1 <- stan("stan/synchrony1.stan", data=c("N","y","J","species","year","type"), iter=1000, chains=4)

library("shinyStan")
launch_shinystan(fit1)

# fit1 <- stan("stan/synchrony1_wpreds.stan", data=c("N","y","J","species","year","type"), iter=1000, chains=4)

# Hack to look at the data
qq<-summary(fit1)
class(qq) # gives the overall summary and then each of the four chains (I think)
qq[[1]][,1]
spptrends.fromstan <- as.vector(qq[[1]][,1])[4:74]
spptrends <- cbind(specieschar.formodel, spptrends.fromstan, lmfits)
spptrends$doybydec <- spptrends$spptrends.fromstan*10
spptrends$doybydec.lm <- spptrends$lmfits*10
spptrends$doybydec.notype <- spptrends.fitnotype.fromstan*10
spptrends$doybydec.wcovar <- spptrends.fitnotypewcovar.fromstan*10

# get a count of years of data for each spp
specieschar.countyrs <- aggregate(rawlong["phenovalue"],
    rawlong[c("intid", "studyid", "species", "spp")], FUN=length)
names(specieschar.countyrs)[names(specieschar.countyrs)=="phenovalue"] <- "nyrs"
specieschar.countyrs.sm <- aggregate(specieschar.countyrs["nyrs"],
    specieschar.countyrs[c("studyid", "species", "spp")], FUN=mean) # lost 20 rows
# merge in nyrs count
spptrends <- merge(spptrends, specieschar.countyrs.sm, by=c("studyid", "species", "spp"))

ggplot(data=spptrends, aes(x=type, y=doybydec))+
geom_path(aes(group=factor(studyid), colour=factor(int_type)))+
    geom_point(aes(colour=factor(int_type), size=1))+ 
    theme_bw()+theme(panel.grid.major = element_blank(),
    axis.title.x =element_text(size=12), axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, angle=90))+
    ylab("Phenological shift (days/decade)")+xlab("Resource (1) or consumer (2)")+
    ggtitle("Slopes from Stan \n resources are on left, consumers on right (I think)")

ggplot(data=spptrends, aes(x=type, y=doybydec.lm))+
geom_path(aes(group=factor(studyid), colour=factor(int_type)))+
    geom_point(aes(colour=factor(int_type), size=1))+ 
    theme_bw()+theme(panel.grid.major = element_blank(),
    axis.title.x =element_text(size=12), axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, angle=90))+
    ylab("Phenological shift (days/decade)")+xlab("Resource (1) or consumer (2)")+
    ggtitle("Slopes from basic lm \n resources are on left, consumers on right (I think)")

# look at effects of type
plot(doybydec~doybydec.lm, spptrends, ylab="Change in days per decade, from Stan",
    xlab="Change in days per decade, simple lm", col="dodgerblue4")
points(doybydec.notype~doybydec.lm, spptrends, pch=1, col="red")
abline(0,1)
abline((lm(doybydec~doybydec.lm, spptrends)), col="dodgerblue4")
abline((lm(doybydec.notype~doybydec.lm, spptrends)), col="red")
legend("topleft", legend=c("with type", "no type", "1:1"), pch=c(1,1,1),
    col=c("dodgerblue4", "red", "black"), bty="n")

# look at effects of wcovar (no type)
plot(doybydec.notype~doybydec.lm, spptrends, ylab="Change in days per decade, from Stan",
    xlab="Change in days per decade, simple lm", col="dodgerblue4")
points(doybydec.wcovar~doybydec.lm, spptrends, pch=1, col="purple")
abline(0,1)
abline((lm(doybydec.notype~doybydec.lm, spptrends)), col="dodgerblue4")
abline((lm(doybydec.notype~doybydec.lm, spptrends)), col="purple")
legend("topleft", legend=c("no type", "no type w covar", "1:1"), pch=c(1,1,1),
    col=c("dodgerblue4", "purple", "black"), bty="n")
posvector <- rep(3, length(spptrends))
text(spptrends$doybydec.lm, spptrends$doybydec.wcovar, labels=as.integer(spptrends$nyrs),
    cex=0.7, pos=posvector)


write.csv(spptrends, "output/synchrony1spptrends.csv", row.names=FALSE)
