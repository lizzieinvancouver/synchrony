// Stole from W. Pearse, here (im Py code) //
// https://github.com/Auerilas/BayesPGLS/blob/master/PGLS.py //

data{
	int N;
	int K;
	vector[N] y;
	matrix[N, N] V;
	matrix[N, N] Lmat;
	vector[N] X;
}
transformed data{
	real Ndiv;
	matrix[N, N] Ident;
	Ndiv <- N;
	Ident <- diag_matrix(rep_vector(1, N));
}
parameters{
	real<lower=0> sigma;
	real<lower=0, upper=1> lambda;
	vector[2] B;
}
transformed parameters{
}
model{
	matrix[N,N] lambda_placeholder;
	matrix[N, N] V_lambda;
	matrix[N, N] V_sig;
	vector[N] yhat;
	real constant;
	real detV;
	real cdf;
	real logLike_PGLS;
	lambda_placeholder <- lambda * Lmat;
	V_lambda <- (lambda_placeholder + Ident) .* V;
	V_sig <- sigma^2*V_lambda;
	yhat <- B[1] + B[2]*X;
	detV <- log_determinant(V_sig);
	cdf <- ((y-yhat)'*inverse_spd(V_sig)*(y-yhat));
	logLike_PGLS <-  -0.5*(detV + cdf);
	increment_log_prob(logLike_PGLS);
	
	B ~ normal(0, 1);
	sigma ~ cauchy(0, 2.5);
	lambda ~ beta(1, 1);
}
