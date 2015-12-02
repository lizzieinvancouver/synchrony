## Started 25 Nov 2015 ##
## By Lizzie ##

## Happy American Thanksgiving! ##

## Random from JD:
# http://www.timetree.org/

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

##
setwd("~/Documents/git/projects/trophsynch/synchrony")

# TOL tree 
# http://datadryad.org/resource/doi:10.5061/dryad.8j60q

library(phytools)

# hinchlifftree1 <- read.newick("input/phylogeny/hinchliff/fromdryad/draftversion3.tre") # too big!

# Dan's approach to the big tree
library(ade4)
phtree <- scan("input/phylogeny/hinchliff/fromdryad/draftversion3.tre", what = "character")
phtre <- newick2phylog(phtree) # Convert it ade4 phylog class
spnames <- names(phtre$leaves) # Get species names

# from supp, see section "Step 4: Prune the Taxonomy to Taxa Required by Phylogenetic Inputs."
hinchlifftree <- read.newick("input/phylogeny/hinchliff/fromSupp/synth-v3-pruned-taxonomy.tre")

mytips <- read.csv("input/specieslist_classfromwiki.csv", header=TRUE)

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

## check for Beth

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

