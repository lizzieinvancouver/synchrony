// Major original edits by M Kosmala
// from file called synchrony1_notype_MK.stan in gelmanhill/synchrony repo
// updated to add covariance matrix by H Kharouba
// 1 December 2015: Commented out random intercepts to see if it improves model fit

data {
  int N;                                // # data points
  int J;                                // # species
  vector[N] y;                          // DOY of pheno event
  int species[N];                       // species identity, coded as int
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
// real mu_a;                            // mean intercept across species
// real<lower=0> sigma_a;                // variation of intercept among species; implicit uniform prior
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
//  a ~ uniform(mu_a, sigma_a); // Tried changing to uniform but "Exception thrown at line // 37: stan::math::uniform_log: Upper bound parameter is 1.21591, but must be greater than // 1.8004" 
  b ~ normal(mu_b, sigma_b);
} 

