data {
  int N;                                  // # data points
  real y[N];
  int J;                                  // # species
  int species[N];
  real year[N];
}
parameters {
 real<lower=0> sigma_y;
  real species_trend[J];                  // trend for species
  real<lower=0> sigma_trend[J];
 real species_level[J];
  real<lower=0> sigma_species_level[J];  
}
model {
  real ypred[N];
  species_level ~ normal(0,1);
  species_trend ~ normal(0,1);
  for (n in 1:N){
    int s;
    s <- species[n];
   ypred[n] <- species_level[s]*sigma_species_level[s]+
                  (sigma_trend[s]*species_trend[s])*year[n];
  }
  y ~ normal(ypred, sigma_y);
} 
