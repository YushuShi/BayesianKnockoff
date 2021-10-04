data {
  int N;           // num obs
  int lengthZ;     // num total covariates
  vector[lengthZ] mu0;                
  matrix[lengthZ,lengthZ] tau0;            
  int Xvec[lengthZ+1,N] ;    
}

parameters {
  vector<lower=-6.0, upper=6.0>[lengthZ] Z; 
  cov_matrix[lengthZ] Sigma; 
}

transformed parameters {
    vector[lengthZ+1] theta;
    vector[lengthZ+1] thetasd;
	theta[1]=1;
	for(i in 1:lengthZ){
		theta[i+1]=exp(Z[i]);
	}
	thetasd=theta/sum(theta);
}

model {
    for(n in 1:N){
		Xvec[,n] ~ multinomial(thetasd);
	}
	Z ~  multi_normal_prec(mu0 , Sigma);
	Sigma ~ wishart(lengthZ,tau0);
}
