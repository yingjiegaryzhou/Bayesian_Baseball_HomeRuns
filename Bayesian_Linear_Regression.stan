data {
  int<lower=0> N;  // number of observ
  vector[N] b_home_run;  // response
  vector[N] player_age;
  vector[N] b_ab;
  vector[N] exit_velocity_avg;
  vector[N] launch_angle_avg;
  int<lower=0> N_new;  // number of new observ for prediction
  vector[N_new] player_age_new;
  vector[N_new] b_ab_new;
  vector[N_new] exit_velocity_avg_new;
  vector[N_new] launch_angle_avg_new;
}
parameters {
  real alpha;  // intercept
  real beta_age;
  real beta_ab;
  real beta_velocity;
  real beta_angle;
  real<lower=0> sigma;  // sd
}
model {
  // Priors
  alpha ~ normal(10, 20);  // weakly informative prior centered around 0
  beta_age ~ normal(0, 1);  // weakly informative prior based on age range
  beta_ab ~ normal(0, 0.1);  // weakly informative prior based on ABs range
  beta_velocity ~ normal(0, 1);  // weakly informative prior based on exit velo range
  beta_angle ~ normal(0, 1);  // weakly informative prior based on launch angle range
  sigma ~ cauchy(0, 5);  // weakly informative prior for sd

  // Likelihood
  b_home_run ~ normal(alpha + beta_age * player_age + beta_ab * b_ab + beta_velocity * exit_velocity_avg + beta_angle * launch_angle_avg, sigma);
}

generated quantities {
  vector[N] b_home_run_sim;  // simulated data
  vector[N_new] b_home_run_pred;  // predictions for new data
  real p_positive_effect_age = step(beta_age);  // prob that age has a positive effect
  real p_positive_effect_ab = step(beta_ab);  // prob that ABs have a positive effect
  real p_positive_effect_velocity = step(beta_velocity);  // prob that exit velo has a positive effect
  real p_positive_effect_angle = step(beta_angle);  // prob that launch angle has a positive effect

  // Generate simulated data and predictions
  for (i in 1:N) {
    b_home_run_sim[i] = normal_rng(alpha + beta_age * player_age[i] + beta_ab * b_ab[i] + beta_velocity * exit_velocity_avg[i] + beta_angle * launch_angle_avg[i], sigma);
  }
  for (i in 1:N_new) {
    b_home_run_pred[i] = normal_rng(alpha + beta_age * player_age_new[i] + beta_ab * b_ab_new[i] + beta_velocity * exit_velocity_avg_new[i] + beta_angle * launch_angle_avg_new[i], sigma);
  }
}

