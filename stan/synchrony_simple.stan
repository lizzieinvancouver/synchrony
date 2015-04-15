// Similar to something once called synchrony1.stan but ...
// removes all reference to trends and fits just an intercept model

data {
  int N;                                  // # data points
  real y[N];
  int J;                                  // # species
  int species[N];
  real year[N];
  int<lower=1,upper=2> type[J];           // species type
}
parameters {
  real<lower=0> sigma_y;
  real level[2];                          // avg level for type
  real species_level[J];
  real<lower=0> sigma_level[2];  
}
model {
  real ypred[N];
  species_level ~ normal(0,1);
  for (n in 1:N){
    int s;
    int t;
    s <- species[n];
    t <- type[s];
    ypred[n] <- level[t] + sigma_level[t]*species_level[s];
  }
  y ~ normal(ypred, sigma_y);
} 
