 //Input Data
data {
  int<lower=2> N; //Num data points							
  vector<lower=0,upper=1>[N] T; //Indicator 0,1 for group
  vector[N] y; //Input log-intensities
}

//Model Parameters
parameters {
  real mu;  //
  real d; //log-e fold-change
  //real tau; //protein variance
  real <lower=0> tau;
}

//Parameters which are transformations of the above
transformed parameters {
  real sigma; //protein standard-deviation (Stan uses sd for normal distr.)
  sigma = sqrt(tau);
}

//Model specification
model {
  mu ~ normal(0,100); //Prior on mu
  d ~ normal(0,100); //Prior on d
  tau ~ inv_gamma(1,0.1); //Prior on tau
  y ~ normal(mu + T*d,sigma); //Likelihood
}


