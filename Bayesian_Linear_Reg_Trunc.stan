// Same as Bayesian_Linear_Regression.stan but with truncated normal likelihood to prevent negative home runs
data {
  int<lower=0> N;
  vector[N] b_home_run;  
  vector[N] player_age;
  vector[N] b_ab;
  vector[N] exit_velocity_avg;
  vector[N] launch_angle_avg;
  int<lower=0> N_new;
  vector[N_new] player_age_new;
  vector[N_new] b_ab_new;
  vector[N_new] exit_velocity_avg_new;
  vector[N_new] launch_angle_avg_new;
}

parameters {
  real alpha;
  real beta_age;
  real beta_ab;
  real beta_velocity;
  real beta_angle;
  real<lower=0> sigma;
}

model {
  // Priors
  alpha ~ normal(10, 20);
  beta_age ~ normal(0, 1);
  beta_ab ~ normal(0, 0.1);
  beta_velocity ~ normal(0, 1);
  beta_angle ~ normal(0, 1);
  sigma ~ cauchy(0, 5);

  // Likelihood
  for (i in 1:N) {
    target += normal_lpdf(b_home_run[i] | alpha + beta_age * player_age[i] + beta_ab * b_ab[i] + beta_velocity * exit_velocity_avg[i] + beta_angle * launch_angle_avg[i], sigma);
    target += -normal_lccdf(0 | alpha + beta_age * player_age[i] + beta_ab * b_ab[i] + beta_velocity * exit_velocity_avg[i] + beta_angle * launch_angle_avg[i], sigma);  // truncation below 0
  } // Truncated normal
}

generated quantities {
  vector[N] b_home_run_sim;
  vector[N_new] b_home_run_pred;
  real p_positive_effect_age = step(beta_age);
  real p_positive_effect_ab = step(beta_ab);
  real p_positive_effect_velocity = step(beta_velocity);
  real p_positive_effect_angle = step(beta_angle);

  // Get simulated data and predictions
  for (i in 1:N) {
    b_home_run_sim[i] = normal_rng(alpha + beta_age * player_age[i] + beta_ab * b_ab[i] + beta_velocity * exit_velocity_avg[i] + beta_angle * launch_angle_avg[i], sigma);
  }
  for (i in 1:N_new) {
    b_home_run_pred[i] = normal_rng(alpha + beta_age * player_age_new[i] + beta_ab * b_ab_new[i] + beta_velocity * exit_velocity_avg_new[i] + beta_angle * launch_angle_avg_new[i], sigma);
  }
}


