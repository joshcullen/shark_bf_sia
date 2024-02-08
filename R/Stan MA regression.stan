
/////////////////////////////////////////////////////////////////
/// Stan code for major axis regression of bite force scaling ///
/////////////////////////////////////////////////////////////////

// Code based on model used in Brose, U., Archambault, P., Barnes, A.D. et al. Predator traits determine food-web architecture across ecosystems. Nat Ecol Evol 3, 919–927 (2019). https://doi.org/10.1038/s41559-019-0899-x


data {
  int n;                     //number of observations
  vector[n] log_FL_c;        //log-transformed, centered fork length
  vector[n] log_ABF;         //log-transformed anterior bite force
}

parameters {
  real log_a;                //log-transformed intercept
  real b;                    //scaling coeff
  real<lower=0> sigma;       //standard deviation of residual errors from Normal distrib
}

model {
  vector[n] xp;
  vector[n] yp;
  vector[n] xd;
  vector[n] yd;
  vector[n] d;

  real epsilon;
  epsilon = 1e-16;           //small non-zero value

  // priors
  log_a ~ normal(0, 5);
  b ~ normal(0, 2);
  sigma ~ gamma(1, 1);

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

  //posterior prediction
  for (i in 1:n){
    mu_pp = log_a + log_FL_c[i] * b;
    y_pp[i] = normal_rng(mu_pp, sigma);
  }

}

