functions {
  real weighted_normal_lpdf(vector y, vector x, real beta1, real beta2, real sigma, vector weights, int n) {
    real lp = 0.0;
    for (i in 1:n) {
      lp += weights[i] * normal_lpdf(y[i] | beta1 + beta2 * x[i], sigma);
    }
    return lp;
  }
  
  real weighted_bernoulli_lpmf(int[] x, real p, vector weights, int n) {
    real lp = 0.0;
    for (i in 1:n) {
      lp += weights[i] * bernoulli_lpmf(x[i] | p);
    }
    return lp;
  }
}

data {
  int<lower=1> n; // sample size
  int<lower=1> n_pred; // predicted sample size  
  vector[n] y;    // Response for linear model
  int x[n];       // Predictor for linear model
  vector<lower=0>[n] weights; // weights for linear model
}

parameters {
  real beta1;
  real beta2;
//  real<lower=0> sigma_beta;
  real<lower=0> sigma;
  real<lower=0> p;
}

//transformed parameters {
  //real mean_x_HT = dot_product(to_vector(x), weights) / total_N; 
  //real lambda = mean_x_HT;
//}

model {
//  beta1 ~ normal(0, 100);  // Priors for linear regression coefficients
//  beta2 ~ normal(0, 100);
//  sigma_beta ~ cauchy(0, 10);
//  sigma ~ cauchy(0, 10);
  
  target += weighted_bernoulli_lpmf(x | p, weights, n);
  target += weighted_normal_lpdf(y | to_vector(x), beta1, beta2, sigma, weights, n);
}

generated quantities {
  vector[n_pred] syn_y;
  int syn_x[n_pred];
  real q_x;
  real q_y;
  for (i in 1:n_pred){
  syn_x[i] = bernoulli_rng(p);
  syn_y[i] = normal_rng(beta1 + beta2 * syn_x[i],sigma);
  }
  q_x = sum(syn_x);
  q_y = sum(syn_y)/n_pred;
}



