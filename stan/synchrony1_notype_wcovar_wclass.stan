// Major original edits by Lizzie
// trying to update synchrony1_notype_wcovar.stan file to add class to model

data {
  int N;                                // # data points
  int J;                                // # species
  int H;				// # class for each species
  vector[N] y;                          // DOY of pheno event
  int species[N];                       // species identity, coded as int
  int class [N];			// class identity
  vector[N] year;                       // year of data point
  int nVars;					//number of predictors
  matrix[nVars, nVars] Imat;
}
parameters {
  vector[J] a;                          // the intercept for each species
  vector[J] b;                          // the slope for each species
  real<lower=0> sigma_y;                // measurement error, noise, etc.
  cov_matrix[nVars] Omega3;
 
  // hyperparameters
  real mu_a;                            // mean intercept across species
  real<lower=0> sigma_a;                // variation of intercept among species; implicit uniform prior
  real mu_b;                            // mean slope across species
  real<lower=0> sigma_b;                // variation of slope among species; implicit uniform prior
}


model {
	real ypred[N];
	//priors for covariance matrix:
	 Omega3 ~ inv_wishart(nVars+1, Imat);  // inverse-Wishart prior for correlations, 	
	 	
  for (i in 1:N){
    ypred[i] <- a[species[i]] + b[species[i]] * year[i];
  }
  y ~ normal(ypred, sigma_y);
  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
} 

