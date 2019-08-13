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
  int<lower=0> M_new;   // number of new observations
  // Omitting random intercepts for the predicted values since the predicted values will be used
  // to compare within subject
  matrix[M_new, K] x_new; // predictor matrix for new observations
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

generated quantities{
  vector[M_new] y_pred; // response vector for new observations
  for(m in 1:M_new){
      y_pred[m] = normal_rng(x_new[m] * beta + alpha, sigma_err);
  }
}
