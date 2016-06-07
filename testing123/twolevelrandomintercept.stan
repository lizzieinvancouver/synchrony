//Two-level (1 hierarchical grouping) random slope model

data{
	int<lower=0> N; 				//Level 1: Number of observations
	int<lower=0> Nspp; 				//Level 2: Number of species (Grouping factor)
	int species[N]; 	//species identity, coded as int
	
	//predictors
	vector[N] year; 	//year of data point
	
	//response
	real y[N]; 		//DOY of pheno event (OR temperature for temperature change model)

}

parameters{	
	// hyperparameters
	real mu_a;				//mean intercept across species (population)
	real mu_b;				//mean slope across species (population)
	real<lower=0> sigma_y; 	//measurement error, noise etc. (population standard deviation)
	
	real a[Nspp]; 		//the intercept for each species (random effect)
	real<lower=0> sigma_a;	//variation of intercept among species; [sd of random effects]

}

transformed parameters{
 //Varying intercepts
 real beta_0[Nspp];
 
 //Individual mean
 real ypred[N];
 
 //Varying intercepts definition
 for (j in 1:Nspp){
 	beta_0[j]<-mu_a+a[j];
 }

//Individual mean
for (i in 1:N){
		ypred[i]<-beta_0[species[i]]+mu_b*year[i];
	}
	
}

model{
	//Random effects distribution	
	a~normal(0, sigma_a);
	
	y~normal(ypred, sigma_y);
}	

	