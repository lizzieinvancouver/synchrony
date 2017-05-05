## Started 5 May 2017 ##
## By Lizzie ##

## Quick plotting by Lizzie ##

library(rgdal)
library(raster)
setwd("~/Documents/git/projects/trophsynch/synchrony_lizzieinvancouver/maps")
medots <- read.csv("~/Documents/git/projects/trophsynch/synchrony_lizzieinvancouver/maps/studies_geog.csv", header=TRUE)
land<-shapefile("input/ne_50m_land/ne_50m_land.shp") 
boundars<-shapefile("input/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp")

## Plotting ...
# plot(land,col="snow1",lty=0)

# color option 1: 
# plot(boundars,col="slategray",border="snow1")
# points(medots$x, medots$y, col = "deeppink", cex = 0.6, lwd=1.5)
# color option 2: 
plot(boundars,col="grey95",border="slategray")
points(medots$x, medots$y, col = "red3", cex = 0.6, lwd=1.5)
