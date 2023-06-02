
/////////////////////////////////////////////////////////////////
/// Stan code for major axis regression of bite force scaling ///
/////////////////////////////////////////////////////////////////


data {
  int n;
  vector[n] log_FL_c;
  vector[n] log_ABF;
}

parameters {
  real log_a;
  real<lower=0> b;
  real<lower=0> sigma;
}

model {
  vector[n] xp;
  vector[n] yp;
  vector[n] xd;
  vector[n] yd;
  vector[n] d;

  real epsilon;
  epsilon = 1e-16;

  // priors
  log_a ~ normal(5, 2.5);
  b ~ normal(2, 1);
  sigma ~ gamma(1,1);

  // likelihood
  // (xp[i],yp[i]) = projection of (x[i],y[i]) on major axis dep. on a,b
  xp = (b * (log_ABF - log_a) + log_FL_c) ./ (b * b + 1);
  yp = log_a + b * xp;

  // residuals
  xd = log_FL_c - xp;
  yd = log_ABF - yp;
  d = sqrt(epsilon + (xd) .* xd + (yd) .* yd);
  0 ~ normal(d, sigma);
}


generated quantities {
  vector[n] y_pp;
  real mu_pp;

  for (i in 1:n){
    mu_pp = log_a + log_FL_c[i] * b;
    y_pp[i] = normal_rng(mu_pp, sigma);
  }

}

