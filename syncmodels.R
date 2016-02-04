### Started 1 April 2015 ###
### April Fool's Day! ###

## Trying to run STAN with synchrony data ##

## Last updated 13 December to organize models a little #
# after looking at interacting species and new model with random intercepts (seems best) ##
## ... and updated 9 Sep 2015 to updated data ##


# Note that Heather calls synchrony1_notype_wcovar.stan 
# synchrony1_notype_hk.stan

rm(list=ls()) 
options(stringsAsFactors = FALSE)

# setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")
setwd("~/Documents/git/projects/trophsynch/synchrony")

library(ggplot2)
library(rstan)
# library(lme4)

# get data
source("datacleaning.R")

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
makehingeplot("177", c("Caterpillar_a spp.", "Parus1 major"))
makehingeplot("178", c("Ficedula1 albicollis"))
makehingeplot("179", c("Quercus2 robur", "Caterpillar_c spp."))
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


stop("Lizzie made the code stop here")


##
## Hinge models
##

# Heather's code to add in data formatted for hinge model
hinge <- subset(rawlong.nodups, intid=="170" | intid=="171" | intid=="177" |
    intid=="178" | intid=="179" | intid=="180" |intid=="181" | intid=="189" |
    intid=="191"|intid=="193" |intid=="194" | intid=="195" |intid=="196"|
    intid=="201" |intid=="207" |intid=="208")
hinge_non <- subset(rawlong.nodups, intid!="170" & intid!="171" & intid!="177"
     & intid!="178" & intid!="179" & intid!="180" & intid!="181" & intid!="189"
     & intid!="191" & intid!="193" & intid!="194" & intid!="195" & intid!="196"
     & intid!="201" & intid!="207" & intid!="208")
hinge_non$newyear<-hinge_non$year
hinge_pre<-subset(hinge, year<=1981); hinge_pre$newyear<-1981
hinge_post<-subset(hinge, year>1981); hinge_post$newyear<-hinge_post$year
hinges<-rbind(hinge_pre, hinge_post)
rawlong.tot<-rbind(hinge_non, hinges)

# prep the data to fit the model including:
# aggregate to get species level characteristics
# subset down to the phenovalues
rawlong.tot$yr1981 <- rawlong.tot$newyear-1981
N <- nrow(rawlong.tot)
y <- rawlong.tot$phenovalue
specieschar.hin<- aggregate(rawlong.tot["phenovalue"], rawlong.tot[c("studyid", "species", "int_type", "spp")], FUN=length)
J <- nrow(specieschar.hin)
species <- as.numeric(as.factor(rawlong.tot$species))
year <- rawlong.tot$yr1981
nVars <-1
Imat <- diag(1, nVars)


##
## New model as of December 2015 here!
## Random slopes only, no type (aka Margaret Kosmala's model) and with hinge
fit.hinge.rs <- stan("stan/synchrony1_notype_randslops_wcovar.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)

library("shinystan")
launch_shinystan(fit.hinge.rs)

# fit1 <- stan("stan/synchrony1_wpreds.stan", data=c("N","y","J","species","year","type"), iter=1000, chains=4)


# Random slopes, random intercepts, no type (aka Margaret Kosmala's model) and with hinge
fit.hinge <- stan("stan/synchrony1_notype_wcovar.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)


###
### New as of December 2015, trying to match up interacting species and look at differences ....
####
fh.sim <- extract(fit.hinge.rs) # extract(fit.notypewcov)
dim(fh.sim$b)
# here's one iteration
fh.sim$b[2000,]

specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species"))

intid <- read.csv("input/raw_aug.csv", header=TRUE)
intid.sm <- subset(intid, select=c("studyid", "spp1", "spp2", "intid" , "interaction"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]

it1000 <- matrix(0, ncol=1000, nrow=192) # should fix someday so NAs go away
## Lots of studyIDs don't match ....
## "HMK006" "HMK007" "HMK010" "HMK012" "HMK015" "HMK020" "HMK021" "HMK027" "HMK031" "HMK048" "KJB003"
for (i in 3000:4000){
    specieschar.formodel.sm$model <- fh.sim$b[i,]
    andtheanswer <- merge(intid.nodups, specieschar.formodel.sm, by.x=c("studyid", "spp1"),
        by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, specieschar.formodel.sm, by.x=c("studyid", "spp2"),
        by.y=c("studyid", "species"), all.x=TRUE)
    it1000[,(i-3000)] <- andtheanswer$model.x-andtheanswer$model.y
}

meanchange <- rowMeans(it1000, na.rm=TRUE)*10

par(mfrow=(c(1,2)))
# get the mean slope shifts
mean(colMeans(fh.sim$b))*10 # they shift about 4 days/decade
hist(colMeans(fh.sim$b)*10, main="", xlab="change in phenology (days/decade)", breaks=20)
hist(meanchange, main="", xlab="change in synchrony (days/decade)", breaks=20)
mean(meanchange, na.rm=TRUE) # they drift apart by half a day a decade

andtheanswer <- cbind(andtheanswer, meanchange)

library(ggplot2)
ggplot(andtheanswer, aes(x=meanchange, fill=interaction)) +
    geom_histogram(binwidth=0.5, alpha=0.5, position="identity")

ggplot(andtheanswer, aes(x=meanchange, colour=interaction)) + geom_density() # I don't like geom_density but my default histogram is awful

# histograms
preds <- subset(andtheanswer, interaction=="predation")
comps <- subset(andtheanswer, interaction=="competition")
polln <- subset(andtheanswer, interaction=="pollination")
herbiv <- subset(andtheanswer, interaction=="herbivory")
paras <- subset(andtheanswer, interaction=="parasitism")
mutul <- subset(andtheanswer, interaction=="mutuliasm")

par(mfrow=c(2,3))
xlim <- c(-17,17)
hist(preds$meanchange, main="predators", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(comps$meanchange, main="competitors", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(polln$meanchange, main="pollination", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(herbiv$meanchange, main="herbivory", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(paras$meanchange, main="parasitism", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(mutul$meanchange, main="mutualism", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)

ggplot(andtheanswer, aes(x=meanchange, fill=interaction)) + 
    geom_histogram(data = , fill = "red", alpha = 0.2) + 
    geom_histogram(data = subset(andtheanswer, interaction==""), fill = "blue", alpha = 0.2) + 
    geom_histogram(data = subset(andtheanswer, interaction==""), fill = "purple", alpha = 0.2) + 
    geom_histogram(data = subset(andtheanswer, interaction=="herbivory"), fill = "green", alpha = 0.2) +
    geom_histogram(data = subset(andtheanswer, interaction=="pollination"), fill = "orange", alpha = 0.2) + 
    geom_histogram(data = subset(andtheanswer, interaction=="parasitism"), fill = "yellow", alpha = 0.2) 

andtheanswer <- cbind(meanchange, andtheanswer)
bigchanges <- subset(andtheanswer, abs(meanchange)>5)
smallerchanges <- subset(andtheanswer, abs(meanchange)<5)

andtheanswer.sm <- subset(andtheanswer, is.nan(meanchange)==FALSE)


##
## Non-hinge models
##

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

# Random slopes, random intercepts, no type model with covar matrix model (aka Margaret Kosmala's model) 
nVars <- 1
Imat <- diag(1, nVars)
fit.notypewcov <- stan("stan/synchrony1_notype_wcovar.stan", data=c("N","y","J","species","year", "nVars", "Imat"), iter=2000, chains=4)
# Also possible to divide the up the model compilation and sampling (per Allen Riddel)
# model <- stan_model("stan/synchrony1_notype_wcovar.stan") 
# fit.notypewcov <- sampling(model, data, iter, chain) # check that sampling is correct command for R

fit.notypewcov <- stan("stan/synchrony1_notype_wcovar.stan", data=c("N","y","J","species","year", "nVars", "Imat"), iter=2000, chains=4)
print(fit.notypewcov)

yy <- summary(fit.notypewcov)
intercepts.fitnotypewcovar.fromstan <- as.vector(yy[[1]][,1])[1:71]
spptrends.fitnotypewcovar.fromstan <- as.vector(yy[[1]][,1])[72:142]


# Random slopes, random intercepts, no type (aka Margaret Kosmala's model)
fit.notype <- stan("stan/synchrony1_notype.stan", data=c("N","y","J","species","year"), iter=2000, chains=3)
print(fit.notype)

zz<-summary(fit.notype)
intercepts.fitnotype.fromstan <- as.vector(zz[[1]][,1])[1:71]
spptrends.fitnotype.fromstan <- as.vector(zz[[1]][,1])[72:142]

# model with type (resource or consumer); the one and only original
fit1 <- stan("stan/synchrony1.stan", data=c("N","y","J","species","year","type"), iter=1000, chains=4)


## Hack to look at the data
## We should check this is CORRECT given updates to Stan and data!
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
