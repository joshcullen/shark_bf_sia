
//////////////////////////////////////////////////////////
/// Stan code for OLS regression of bite force scaling ///
//////////////////////////////////////////////////////////


data {
  int n;                     //number of observations
  vector[n] log_FL_c;        //log-transformed, centered fork length
  vector[n] ABF;             //anterior bite force
}

parameters {
  real log_a;                //log-transformed intercept
  real b;                    //scaling coeff
  real<lower=0> sigma;       //standard deviation of residual errors from Normal distrib
}

model {
  vector[n] log_mu;
  vector[n] mu;

  // priors
  log_a ~ normal(5, 1);
  b ~ normal(2, 1);
  sigma ~ normal(0, 1);

  // likelihood
  log_mu = log_a + log_FL_c * b;
  mu = exp(log_mu);
  ABF ~ normal(mu, sigma);
}


generated quantities {
  vector[n] y_pp;
  real log_mu_pp;
  real mu_pp;

  //posterior prediction
  for (i in 1:n){
    log_mu_pp = log_a + log_FL_c[i] * b;
    mu_pp = exp(log_mu_pp);
    y_pp[i] = normal_rng(mu_pp, sigma);
  }

}
