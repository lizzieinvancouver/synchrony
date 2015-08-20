// Code from meeting with Andrew 23 Apr 2015
// But updated by Lizzie to try to get simple model to run

data {
  int N;                                  // # data points
  real y[N];
  int J;                                  // # species
  int species[N];
  real year[N];
  int<lower=1,upper=2> type[J];           // species type
}
parameters {
  real trend[2];                          // avg trend for type
  real<lower=0> sigma_y;
  real species_trend[J];                  // trend for species, relative to trend for type
  real<lower=0> sigma_trend[2];
  real level[2];                          // avg level for type
  real species_level[J];
  real<lower=0> sigma_level[2]; 
  vector[J] a;
  vector[J] b; 
  real<lower=0> sigma_y_dumb;
}
model {
  real ypred[N];
  real ypred_dumb[N];
  species_level ~ normal(0,1);
  species_trend ~ normal(0,1);
  for (n in 1:N){
    int s;
    int t;
    s <- species[n];
    t <- type[s];
    ypred[n] <- level[t] + sigma_level[t]*species_level[s] +
                  (trend[t] + sigma_trend[t]*species_trend[s])*year[n];
    ypred_dumb[n] < a[s] + b[s]*year[n];
  }
  a ~ normal(0, 1e3);
  b ~ normal(0, 1e3);
  y ~ normal(ypred, sigma_y);
  y ~ normal(ypred_dumb, sigma_y_dumb);
} 
