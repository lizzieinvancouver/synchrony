## Started 25 Nov 2015 ##
## By Lizzie ##

## Happy American Thanksgiving! ##

## Trying to look at effects of phylogeny for synchrony work ##
# Needs one last task done, then I think we can say phylogeny would not change much #
# And thus it is probably not worth the effort (which would be a lot I think) #
# to add in phylogeny ##

## Random from JD: A way (not so good, but something) to calculate distances: #
# http://www.timetree.org/

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

##
setwd("~/Documents/git/projects/trophsynch/synchrony")

## get the species info
# This file could use a little work -- I did it quickly so I did not probably pick the best #
# closely-related taxa (e.g., would be good to pick freshwater Daphnia species when needed etc. #
# I pulled OTTs from the synonyms file (just easier, but could use taxonomy file) #
# Also, I just noticed Cyclops vicinus I did wrong ##

mytips <- read.csv("input/specieslist_classfromwiki.csv", header=TRUE)

##
## Progress from 11 December 2015 (thanks to JD)
## see also: https://cran.r-project.org/web/packages/rotl/vignettes/how-to-use-rotl.html

## update from late Jan 2016
# Yes, rotl seems broken, try installing from github
# https://github.com/ropensci/rotl
# This seems to fix things (I had to specify build_vignette=F to get it to install)
# install_github("ropensci/rotl", dependencies = TRUE, build_vignette=FALSE)
# Note that install_github requires devtools package
# As of 14 April 2016 I got the below to run, but I had to update R first!
library(rotl)
library(ape)
tr <- tol_induced_subtree(ott_ids=mytips$ott_id)
tr <- compute.brlen(tr, method = "Grafen") # a very poor way to calculate branch lengths, we should impove it if we think phylogeny might matter at all!

dat<-read.csv("output/synchrony1spptrends.csv")
lmfits <- read.csv("input/lmfits.csv", header=TRUE) # new file as of April 2016 with error from fits

ott <- tr$tip.label
pat <- "^.*?_ott([0-9]*)"
test <- sub(pat, "\\1", ott)
my.sp <- mytips[match(test, as.character(mytips$ott_id)),]

match.tree<-tr
match.tree$tip.label<-as.vector(my.sp$species)

# look for phylogentic signal in slopes, it would be much better to do this with the variance
library(picante)
 k.dat<-dat$lmfits
 names(k.dat)<-dat$species
phylosignal(k.dat, multi2di(match.tree), reps = 999)

 k.var<-lmfits$varfits
 names(k.var)<-lmfits$uniquespp
phylosignal(k.var, multi2di(match.tree), reps = 999)

# pretty plots
library(phytools)
library(picante)

plot.slopes<-match.phylo.data(match.tree, k.dat)
contMap(plot.slopes$phy, plot.slopes$data)


## TOL trees ##
# Pulled from the Hinchliff paper, I never got the large one (draftversion3.tre) #
# via any method, even using a big server style computer, J Davies found the above method, which #
# is better ##
# http://datadryad.org/resource/doi:10.5061/dryad.8j60q

# hinchlifftree1 <- read.newick("input/phylogeny/hinchliff/fromdryad/draftversion3.tre") # too big!

# Dan's approach to the big tree
library(ade4)
phtree <- scan("input/phylogeny/hinchliff/fromdryad/draftversion3.tre", what = "character")
phtre <- newick2phylog(phtree) # Convert it ade4 phylog class
spnames <- names(phtre$leaves) # Get species names

# from supp, see section "Step 4: Prune the Taxonomy to Taxa Required by Phylogenetic Inputs."
hinchlifftree <- read.newick("input/phylogeny/hinchliff/fromSupp/synth-v3-pruned-taxonomy.tre")

tree<-drop.tip(hinchlifftree, hinchlifftree$tip.label[!hinchlifftree$tip.label %in% paste("ott", mytips$ott_id, sep="")])

## change in tip names and plot
tipnames<-tree$tip.label
mytips$ott_idwott <- paste("ott", mytips$ott_id, sep="")
getnames <- mytips[which(mytips$ott_idwott %in% tipnames),]

getnamesnodup <- getnames[!duplicated(getnames$ott_idwott),]
matchy <- getnamesnodup[match(tipnames, getnamesnodup$ott_idwott),]

tree2<-tree
tree2$tip.label<-as.vector(matchy$species)
plot(tree2)

## check which tip labels are not in tree

for (i in seq_along(unique(mytips$ott_idwott))){
    print(i)
    print(grep(unique(mytips$ott_idwott)[i], phtree))
}

# not found (4 total) on first pass, now fixed
grep("ott4311804", phtree)
grep("ott5521218", phtree)
grep("ott4921563", phtree)
grep("ott4916717", phtree)

## Make a VCV
library(geiger)
tree2$tip.label <- as.vector(matchy$species)
tree2$edge.length <- rep(1,length(tree$edge[,1]))
unclass(tree2)
plot(tree2)
my.vcv <- vcv(tree2)

# https://stat.ethz.ch/pipermail/r-sig-phylo/2014-March/003364.html

