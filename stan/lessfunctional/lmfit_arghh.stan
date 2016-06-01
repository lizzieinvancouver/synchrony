data { 
    int M; 
    int N[M]; 
    matrix[M,N] y; 
  } 
  parameters { 
    real mu[M]; 
  } 
  model { 
    for (m in 1:M) 
      y[m] ~ normal(mu[m], 1); 
  }
