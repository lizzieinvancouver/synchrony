### Started 17 March 2015 ###
### St. Patrick's Day! ###

## Plots with Heather's synchrony data ##
## Some bits adapted from Heather's code ##

options(stringsAsFactors = FALSE)

setwd("~/Documents/git/projects/trophsynch/synchrony")
library(ggplot2)
library(reshape)

# in the data file: spp1 = neg= negative species e.g. resource 
data <- read.csv("input/alldata_lizzie.csv", header=TRUE)
studies <- read.csv("input/studies.csv", header=TRUE)

studiessm <- subset(studies, select=c("studyid", "terrestrial2"))

# clean raw data
source("source/additionalcleaning_lizzie.R")
raw <- mat
write.csv(raw, "output/raw.csv", row.names=FALSE)

# reshape into long format (do a double melt)
raw.sm <- subset(raw, select=c("intid", "year", "studyid", "spp1", "spp2",
    "interaction", "taxatype", "neg_phenovalue", "neg_slope",
    "pos_phenovalue", "pos_slope"))
raw.melt <- melt(raw.sm, id.var=c("intid", "year", "studyid", "interaction",
    "taxatype", "spp1", "spp2"), measure.var=c("neg_phenovalue",
    "pos_phenovalue", "pos_slope", "neg_slope"))

names(raw.melt)[names(raw.melt)=="value"] <- "measurement"
names(raw.melt)[names(raw.melt)=="variable"] <- "whatmeasured"
raw.mm <- melt(raw.melt, id.var=c("intid", "year", "studyid", "interaction",
    "taxatype", "whatmeasured", "measurement"), measure.var=c("spp1",
    "spp2"))
# relabel the pos and neg since they are no longer needed
raw.mm$whatmeasured <- as.character(raw.mm$whatmeasured)
raw.mm[grep("slope", raw.mm$whatmeasured),]$whatmeasured <- "slope"
raw.mm[grep("phenovalue", raw.mm$whatmeasured),]$whatmeasured <- "phenovalue"

write.csv(raw.mm, "output/alldata_long.csv", row.names=FALSE)


# I used a loop to calculate spp pheno shifts only over the years where both species had pheno data.
# I want the interaction type and the minyear cols also
# Need to also add colors by biome or habitat

rowcount<-1
datawmax <- aggregate(data["year"], by=data[c("studyid","intid")], FUN=max); names(datawmax)[3]<-"maxyear"
datawmax <- merge(data, datawmax, by=c("studyid","intid"), all.x=TRUE)
datawmaxwbiome <- merge(datawmax, studiessm, by="studyid", all.x=TRUE)
mat<-subset(datawmaxwbiome, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
Bgroups<-unique(mat$intid); b<-Bgroups; b<-as.character(b)
nup<-data.frame(array(0, c(80, 8)))
names(nup)[1] <- "studyid"; names(nup)[2]<-"species"; names(nup)[3]<-"intid";
    names(nup)[4] <- "interaction"; names(nup)[5] <- "minyear";  names(nup)[6] <- "maxyear";
    names(nup)[7] <- "biome"; names(nup)[8]<-"doybyyr"
mode(nup$interaction) <- "character"
for(i in 1:length(b)) {
yo<-mat[mat$intid==b[i],]
nup[rowcount,1]<-yo[1,c("studyid")]
nup[rowcount+1,1]<-yo[1,c("studyid")]
nup[rowcount,2]<-yo[1,c("spp1")]
nup[rowcount+1,2]<-yo[1,c("spp2")]
nup[rowcount,3]<-yo[1,c("intid")]
nup[rowcount+1,3]<-yo[1,c("intid")]
nup[rowcount,4]<-as.character(yo[1,c("interaction")])
nup[rowcount+1,4]<-as.character(yo[1,c("interaction")])
nup[rowcount,5]<-yo[1,c("minyear")]
nup[rowcount+1,5]<-yo[1,c("minyear")]
nup[rowcount,6]<-yo[1,c("maxyear")]
nup[rowcount+1,6]<-yo[1,c("maxyear")]
nup[rowcount,7]<-as.character(yo[1,c("terrestrial2")])
nup[rowcount+1,7]<-as.character(yo[1,c("terrestrial2")])
m1<-with(yo, lm(neg_phenovalue~year))
m2<-with(yo, lm(pos_phenovalue~year))
nup[rowcount,8]<-summary(m1)$coefficients[2]
nup[rowcount+1,8]<-summary(m2)$coefficients[2]
rowcount<-rowcount+2
}
nup$spp<-with(nup, c(1:2)) 
unique(nup$interaction)



###############
#### Graphs ###
###############

# Plot the raw time series data
pdf("graphs/realdata.pdf", height=4, width=6)
colors=c("blue", "red")
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(neg_phenovalue~year, data=raw, type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw real data (resources and consumers are blue and red)")
for (j in 1:length(unique(raw$spp1))){
  subbydata <- subset(raw, spp1==unique(raw$spp1)[j])
  lines(neg_phenovalue~year, data=subbydata, col=colors[1])
}
for (j in 1:length(unique(raw$spp2))){
  subbydata <- subset(raw, spp2==unique(raw$spp2)[j])
  lines(pos_phenovalue~year, data=subbydata, col=colors[2])
}
dev.off()



# idea to trick plot into what we want:
# replace spp 1 OR 2 in the spp col with the minyear
nup$doybydec <- nup$doybyyr*10
nup$spphack <- ifelse(nup$spp<2, nup$minyear, nup$maxyear)
# plot with interaction types
ggplot(data=nup, aes(x=spphack, y=doybydec))+
geom_path(aes(group=factor(intid), colour=factor(interaction)), size=1)+
    theme_bw()+theme(panel.grid.major = element_blank(),
    axis.title.x =element_text(size=12), axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, angle=90))+
    ylab("Phenological shift (days/decade)")+xlab("Start to end year of study")+
    ggtitle("48 interactions across 25 studies\nresources are on left, consumers on right")

# plot with biomes
ggplot(data=nup, aes(x=spphack, y=doybydec))+
geom_path(aes(group=factor(intid), colour=factor(biome)), size=1)+
    theme_bw()+theme(panel.grid.major = element_blank(),
    axis.title.x =element_text(size=12), axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, angle=90))+
    ylab("Phenological shift (days/decade)")+xlab("Start to end year of study")+
    ggtitle("48 interactions across 25 studies\nresources are on left, consumers on right")

# simple plot
ggplot(data=nup, aes(x=spp, y=doybydec))+
geom_path(aes(group=factor(intid), colour=factor(interaction)), size=1)+
    geom_point(aes(colour=factor(interaction), size=1))+ 
    theme_bw()+theme(panel.grid.major = element_blank(),
    axis.title.x =element_text(size=12), axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, angle=90))+
    ylab("Phenological shift (days/decade)")+xlab("Resource (1) or consumer (2)")+
    ggtitle("48 interactions across 25 studies\nresources are on left, consumers on right")

# simple plot with biome
ggplot(data=nup, aes(x=spp, y=doybydec))+
geom_path(aes(group=factor(intid), colour=factor(biome)), size=1)+
    geom_point(aes(colour=factor(biome), size=1))+ 
    theme_bw()+theme(panel.grid.major = element_blank(),
    axis.title.x =element_text(size=12), axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, angle=90))+
    ylab("Phenological shift (days/decade)")+xlab("Resource (1) or consumer (2)")+
    ggtitle("48 interactions across 25 studies\nresources are on left, consumers on right")
