### Started 4 May 2014 ###
### By Lizzie ###

## For Heather Kharouba's synchrony project ##
## Adapted from timeseries.R (tempeco project) ##


##
## make random data 
##
shabam <- function(curvenum, output, errormultiplier, trend) {
curvey <- seq(from=-(curvenum)*pi, to=curvenum*pi, length.out=output)
trendy <- seq(from=0, to=1, length.out=output)*trend
dater <- sin(curvey)+trendy+(errormultiplier*rnorm(output))
return(dater) }

addevents <- function(dater, eventsnow, magnitudes) {
eventsy <- rep(0, length(dater))
eventsy[eventsnow] <- magnitudes
daterout <- dater+eventsy
return(daterout)
}

# example
goober <- shabam(5, 200, 0.35, 10)
goo <- addevents(goober, c(10,100,50), rep(10,5))
plot(goo, type="l")

# for figure in temporal ecology ms
# 3-part one #
oldendays <- shabam(6, 600, 0.35, 0)
nowadays <- addevents(oldendays, c(20, 130, 220, 330, 440, 570),
    c(1.6, 2.8, 2.9, 2.6, 2, 3.2, -2))
futureshop <- shabam(6, 600, 0.35, 7)
futureshop.e <- addevents(futureshop, c(20, 130, 220, 330, 440, 570),
    c(1.6, 2.8, 2.9, 2.6, 2, 3.2, -2))

par(ps=9)
quartz(width=3.5, height=6)
par(mfrow=c(3,1))
par(oma=c(2,2,2,1)) # set margins 
par(mar=c(2,2,1,2)) # set margins
plot(oldendays, type="l", ylim=c(-3,9))
plot(nowadays, type="l", ylim=c(-3,9))
plot(futureshop.e, type="l", ylim=c(-3,9))


## make data with nonstationarity #
## after a certain period ##

stationary <- shabam(4, 700, 0.25, 0)
nonstat <- shabam(3, 600, 0.25, 4)

stationary1 <- shabam(4, 720, 0.15, 0)
nonstat1 <- shabam(3, 580, 0.15, 4)

stationaryalways <- shabam(8, 1400, 0.15, 0)

quartz(width=6, height=4)
par(mfrow=c(1,1))
nonstat <- c(stationary, nonstat)
plot(nonstat, type="l")
nonstat1 <- c(rep(0, 30), stationary1, nonstat1)
# lines(nonstat1, col="red")
lines(stationaryalways, col="blue")


