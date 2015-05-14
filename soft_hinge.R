soft_hinge <- function(x,a,b){
  z <- (x-b)/a
  return(z/(1-exp(-z)))
}
curve(soft_hinge(x, a=1, b=1975), from=1961, to=2010)

