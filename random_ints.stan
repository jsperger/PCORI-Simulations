data {
  int<lower=0> M;   // number of observations
  int<lower=0> N; // number of subjects
  //int<lower=1> T; // number of treatments
  
  int<lower=0> K;   // number of predictors (including treatments, not including the intercept)
  row_vector[K] x[M];   // predictor matrix
  int<lower=1> subj[M]; // subject id
  vector[M] y;      // outcome vector
  
  // Priors
  real a; // prior intercept
  real<lower=0> a_var; //intercept prior variance
  vector[K] b; // beta prior
  vector<lower=0>[K] b_var; // beta prior variance
  
  // For predictions
  //int<lower=0> M_new;   // number of observations
  //int<lower=0> N_new; // number of subjects
  //matrix[M_new, K] x_new; // predictor matrix for new observations
  //vector[M_new] y_new; // response vector for new observations
  
}

parameters {
  real alpha;           // fixed intercept
  vector[N] b_i;        // random intercepts
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma_err;  // error sd
  real<lower=0> sigma_b; // between subject sd
}
model {
  // priors
  alpha ~ normal(a, a_var);
  beta ~ normal(b, b_var);
  b_i ~ normal(0, sigma_b);
  
  // Likelihood
  for(m in 1:M){
    y[m] ~ normal(x[m] * beta + alpha + b_i[subj[m]], sigma_err);
  }
}

