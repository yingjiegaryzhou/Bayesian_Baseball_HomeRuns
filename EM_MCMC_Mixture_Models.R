##Stan Two Component Mixture

#Load the necessary libraries
library(Lahman)
library(rstan)
library(rstanarm)
library(tidyverse)
library(bayesplot)
library(ggplot2)
library(Lahman)
data("Batting")

#### Test with manual MCMC ####
# Subset the Batting data for seasons after 1975
season_all_relevant <- subset(Batting, yearID > 1975)

allseason_min_atbat <- season_all_relevant %>% filter(AB >= 100)
allseason_min_ab <- allseason_min_atbat %>% filter(HR>0)
hr <- allseason_min_ab$HR
ab <- allseason_min_ab$AB
hrprop <- hr/ab
Y <- hrprop
n <- length(Y)
year <- allseason_min_ab$yearID
player <- allseason_min_ab$playerID

# starting point
alpha <- 0.1
mu0 <- 0.02
mu1 <- 0.06
sigsq0 <- 0.005
sigsq1 <- 0.01

# gibbs sampler
numiters <- 2000
params <- matrix(NA, nrow=numiters, ncol=5)
params[1,] <- c(alpha, mu0, mu1, sigsq0, sigsq1)
inds <- matrix(NA, nrow=numiters-1, ncol=n)
for (i in 2:numiters){
  ## sampling indicator variables
  I <- rep(NA,n)
  for (j in 1:n){
    a <- dnorm(Y[j], mu1, sqrt(sigsq1))*alpha    
    b <- dnorm(Y[j], mu0, sqrt(sigsq0))*(1-alpha)
    p <- a/(a+b)
    I[j] <- rbinom(1, 1, p)
  }
  ## calculating statistics from indicator variables
  n0 <- sum(I==0)
  n1 <- sum(I==1)
  meanY0 <- mean(Y[I==0])
  meanY1 <- mean(Y[I==1])
  ## sampling alphas
  alpha <- rbeta(1,n1+1,n0+1)
  ## sampling means
  mu0 <- rnorm(1,meanY0,sqrt(sigsq0/n0))
  mu1 <- rnorm(1,meanY1,sqrt(sigsq1/n1))
  ## calculating statistics from new means 
  ss0 <- sum((Y[I==0]-mu0)^2)
  ss1 <- sum((Y[I==1]-mu1)^2)
  ## sampling variances
  temp <- rgamma(1, shape=(n0/2)+1, rate=ss0/2)
  sigsq0 <- 1/temp
  temp <- rgamma(1, shape=(n1/2)+1, rate=ss1/2)
  sigsq1 <- 1/temp
  ## storing current values
  params[i,] <- c(alpha, mu0, mu1, sigsq0, sigsq1)
  inds[i-1,] <- I
  print(i)
  #print(round(c(alpha,mu0,mu1,sigsq0,sigsq1),5))
}

params.data1 <- params 
params.data1

#running Gibbs sampler for alternate starting point with code above
alpha <- 0.80
mu0 <- 0.03
mu1 <- 0.05
sigsq0 <- 0.005
sigsq1 <- 0.003

for (i in 2:numiters){
  ## sampling indicator variables
  I <- rep(NA,n)
  for (j in 1:n){
    a <- dnorm(Y[j], mu1, sqrt(sigsq1))*alpha    
    b <- dnorm(Y[j], mu0, sqrt(sigsq0))*(1-alpha)
    p <- a/(a+b)
    I[j] <- rbinom(1, 1, p)
  }
  ## calculating statistics from indicator variables
  n0 <- sum(I==0)
  n1 <- sum(I==1)
  meanY0 <- mean(Y[I==0])
  meanY1 <- mean(Y[I==1])
  ## sampling alphas
  alpha <- rbeta(1,n1+1,n0+1)
  ## sampling means
  mu0 <- rnorm(1,meanY0,sqrt(sigsq0/n0))
  mu1 <- rnorm(1,meanY1,sqrt(sigsq1/n1))
  ## calculating statistics from new means 
  ss0 <- sum((Y[I==0]-mu0)^2)
  ss1 <- sum((Y[I==1]-mu1)^2)
  ## sampling variances
  temp <- rgamma(1, shape=(n0/2)+1, rate=ss0/2)
  sigsq0 <- 1/temp
  temp <- rgamma(1, shape=(n1/2)+1, rate=ss1/2)
  sigsq1 <- 1/temp
  ## storing current values
  params[i,] <- c(alpha, mu0, mu1, sigsq0, sigsq1)
  inds[i-1,] <- I
  print(i)
  #print(round(c(alpha,mu0,mu1,sigsq0,sigsq1),5))
}

params.data2 <- params 
params.data2

# comparing iterations between two chains 
par(mfrow=c(1,1))
plot(1:numiters, params.data2[,1], main="alpha", type="l", col=3, xlab="Number of Iterations", ylab="alpha")
lines(1:numiters,params.data1[,1],col=2)

par(mfrow=c(3,2))
ymin <- min(params.data2[,1],params.data1[,1])
ymax <- max(params.data2[,1],params.data1[,1])
plot(1:numiters,params.data2[,1],main="alpha",type="l",col=3,ylim=c(ymin,ymax), xlab="Number of Iterations", ylab="alpha")
lines(1:numiters,params.data1[,1],col=2)

ymin <- min(params.data2[,2],params.data1[,2])
ymax <- max(params.data2[,2],params.data1[,2])
plot(1:numiters,params.data2[,2],main="mu1",type="l",col="green",ylim=c(ymin,ymax), xlab="Number of Iterations", ylab="mu1")
lines(1:numiters,params.data1[,2],col=2)

ymin <- min(params.data2[,3],params.data1[,3])
ymax <- max(params.data2[,3],params.data1[,3])
plot(1:numiters,params.data2[,3],main="mu2",type="l",col="green",ylim=c(ymin,ymax), xlab="Number of Iterations", ylab="mu2")
lines(1:numiters,params.data1[,3],col=2)

ymin <- min(params.data2[,4],params.data1[,4])
ymax <- max(params.data2[,4],params.data1[,4])
plot(1:numiters,params.data2[,4],main="sigsq1",type="l",col="green",ylim=c(ymin,ymax), xlab="Number of Iterations", ylab="sigsq1")
lines(1:numiters,params.data1[,4],col=2)

ymin <- min(params.data2[,5],params.data1[,5])
ymax <- max(params.data2[,5],params.data1[,5])
plot(1:numiters,params.data2[,5],main="sigsq2",type="l",col="green",ylim=c(ymin,ymax), xlab="Number of Iterations", ylab="sigsq2")
lines(1:numiters,params.data1[,5],col=2)

params.data1.postburn<-params.data1[201:2000,]
params.data2.postburn<-params.data2[201:2000,]

par(mfrow=c(3,2))
acf(params.data1.postburn[,1],lag.max=100,main="ACF: alpha, chain 1")
acf(params.data2.postburn[,1],lag.max=100,main="ACF: alpha, chain 2")
acf(params.data1.postburn[,2],lag.max=100,main="ACF: mu1, chain 1")
acf(params.data2.postburn[,2],lag.max=100,main="ACF: mu1, chain 2")
acf(params.data1.postburn[,3],lag.max=100,main="ACF: mu2, chain 1")
acf(params.data2.postburn[,3],lag.max=100,main="ACF: mu2, chain 2")

par(mfrow=c(2,2))
acf(params.data1.postburn[,4],lag.max=100,main="ACF: sigsq1, chain 1")
acf(params.data2.postburn[,4],lag.max=100,main="ACF: sigsq1, chain 2")
acf(params.data1.postburn[,5],lag.max=100,main="ACF: sigsq2, chain 1")
acf(params.data2.postburn[,5],lag.max=100,main="ACF: sigsq2, chain 2")

# taking only every fiftieth draw 
temp <- 50*c(1:(1800/50))

params.data1.thinned <- params.data1.postburn[temp,]
params.data2.thinned <- params.data2.postburn[temp,]

par(mfrow=c(3,2))
acf(params.data1.thinned[,1],lag.max=100,main="ACF: alpha, chain 1")
acf(params.data2.thinned[,1],lag.max=100,main="ACF: alpha, chain 2")
acf(params.data1.thinned[,2],lag.max=100,main="ACF: mu1, chain 1")
acf(params.data2.thinned[,2],lag.max=100,main="ACF: mu1, chain 2")
acf(params.data1.thinned[,3],lag.max=100,main="ACF: mu2, chain 1")
acf(params.data2.thinned[,3],lag.max=100,main="ACF: mu2, chain 2")

par(mfrow=c(2,2))
acf(params.data1.thinned[,4],lag.max=100,main="ACF: sigsq1, chain 1")
acf(params.data2.thinned[,4],lag.max=100,main="ACF: sigsq1, chain 2")
acf(params.data1.thinned[,5],lag.max=100,main="ACF: sigsq2, chain 1")
acf(params.data2.thinned[,5],lag.max=100,main="ACF: sigsq2, chain 2")

#Combining chains and calculating posterior intervals
params.final <- rbind(params.data1.thinned,params.data2.thinned)

numsamples <- length(params.final[,1])
alpha.gibbs <- params.final[,1]
mu0.gibbs <- params.final[,2]
mu1.gibbs <- params.final[,3]
sigsq0.gibbs <- params.final[,4]
sigsq1.gibbs <- params.final[,5]

par(mfrow=c(3,2))
hist(mu0.gibbs, main="mu1")
hist(mu1.gibbs, main="mu2")
hist(sigsq0.gibbs, main="sigsq1")
hist(sigsq1.gibbs, main="sigsq2")
hist(alpha.gibbs, main="alpha")

summary(mu0.gibbs)
summary(mu1.gibbs)
summary(sigsq0.gibbs)
summary(sigsq1.gibbs)
summary(alpha.gibbs)

#############################################
#Plot equal-variances fitted mixture density
alpha <- finalparam[1]
mu1 <- finalparam[2]
mu2 <- finalparam[3]
sigsq <- finalparam[4] 
par(mfrow=c(1,1))
hist(hrprop, prob=T, ylim = c(0,25))
x <- ppoints(1000)*0.15
y1 <- (1-alpha)*dnorm(x,mu1,sqrt(sigsq))
y2 <- alpha*dnorm(x,mu2,sqrt(sigsq))
y<- (1-alpha)*dnorm(x,mu1,sqrt(sigsq)) + alpha*dnorm(x,mu2,sqrt(sigsq)) 
lines(x,y1,col='green')
lines(x,y2,col='red')
lines(x,y,col='purple', lwd=2,lty=2)
legend('topright', legend=c('Non-elite group', 'Elite group', 'Mixture'),
        col=c('green','red', 'purple'), lty=c(1,1,2), lwd=c(1,1,2))
#############################################

#Individual HR Probabilities for Players
prob.elite.hrprop <- Estep(hrprop, alpha, mu1, mu2, sigsq)

#Find out the playerIDs of player seasons that have a greater than 50% chance of elite
eliteplayernames <- playerID[which(prob.elite.hrprop > 0.50)]

#Find out the years of player seasons that have a greater than 50% chance of elite
eliteplayeryears <- yearID[which(prob.elite.hrprop > 0.50)]

#Find out the home run rates and totals of player seasons that have a elite seasons
eliteplayerhrprop <- hrprop[which(prob.elite.hrprop > 0.50)]
eliteplayerhrtot <- HR[which(mu1.gibbs > 0.50)]
length(mu1.gibbs)

#Get the probability of being in the elite group for player seasons
prob.elite <- Estep(eliteplayerhrprop, alpha, mu1, mu2, sigsq)

#Combine columns and put into table
elite.hrhitters <- cbind(eliteplayernames, eliteplayeryears, prob.elite, eliteplayerhrprop)
final_table <- xtable(elite.hrhitters)
names(final_table) <- c("Player Names", "Year", "Elite Group Probability", "HR Prop")

#Filter out those players who have had seasons of elite home run hitting 
final_table %>% arrange(desc(prob.elite)) %>% count()
final_table %>% arrange(desc(prob.elite)) %>% head(5)
final_table %>% arrange(prob.elite) %>% head(5)

################Stan Implementation################
# Subset the Batting data for seasons after 1975
season_all_relevant <- subset(Batting, yearID > 1975)
allseason_min_atbat <- season_all_relevant %>% filter(AB >= 100)
allseason_min_ab <- allseason_min_atbat %>% filter(HR>0)
hr <- allseason_min_ab$HR
ab <- allseason_min_ab$AB
hrprop <- hr/ab
Y <- hrprop
n <- length(Y)
year <- allseason_min_ab$yearID
player <- allseason_min_ab$playerID

# Define the Stan model
stan_code <- "
data {
    int<lower=0> n; // Number of data points
    vector[n] Y; // Data
}
parameters {
    real<lower=0, upper=1> alpha; // Mixing proportion
    ordered[2] mu; // Means of the two components
    real<lower=0> sigma[2]; // Standard deviations of the two components
}
model {
    alpha ~ beta(1, 1);
    mu[1] ~ normal(0.06, 1);
    mu[2] ~ normal(0.01, 1);
    sigma ~ inv_gamma(2, 0.1);
    for (i in 1:n) {
        target += log_mix(alpha,
                          normal_lpdf(Y[i] | mu[1], sigma[1]),
                          normal_lpdf(Y[i] | mu[2], sigma[2]));
    }
}
"

#Add generated quantities
stan_code <- "
data {
    int<lower=0> n; // Number of data points
    vector[n] Y; // Data
}
parameters {
    real<lower=0, upper=1> alpha; // Mixing proportion
    ordered[2] mu; // Means of the two components
    real<lower=0> sigma[2]; // Standard deviations of the two components
}
model {
    alpha ~ beta(1, 1);
    mu[1] ~ normal(0.06, 1);
    mu[2] ~ normal(0.01, 1);
    sigma ~ inv_gamma(2, 0.1);
    for (i in 1:n) {
        target += log_mix(alpha,
                          normal_lpdf(Y[i] | mu[1], sigma[1]),
                          normal_lpdf(Y[i] | mu[2], sigma[2]));
    }
}
generated quantities {
  vector[n] Y_hat;
  for (i in 1:n) {
    Y_hat[i] = normal_rng(mu[1], sigma[1]) * alpha + normal_rng(mu[2], sigma[2]) * (1 - alpha);
  }
}
"

options(mc.cores = parallel::detectCores())
rstan::rstan_options(mc.cores = parallel::detectCores())

#Compile
stan_model <- stan_model(model_code = stan_code)

#Fit model
fit <- sampling(stan_model, data = list(n = n, Y = Y), iter = 3000, chains = 4)

setwd("/Users/garyzhou/Downloads")
#save(fit, file = "trunc_norm_mixture_mod.rda")
#save(fit, file = "mixture_mod_fit_pp.rda")
#save(fit, file = "mixture_mod_fit.rda")

# Print the fit
print(fit)

# Posterior predictive checks
pp_check(fit)

# Trace plots
mcmc_trace(fit, pars = c("alpha", "mu[1]", "mu[2]", "sigma[1]", "sigma[2]"))

# Autocorrelation plots
mcmc_acf_bar(fit, pars = c("alpha", "mu[1]", "mu[2]", "sigma[1]", "sigma[2]"))
mcmc_acf(fit, pars = c("alpha", "mu[1]", "mu[2]", "sigma[1]", "sigma[2]"))

##Generate posterior predictive samples
#PP Check Test
posterior_samples_test <- as.matrix(fit)
simulated_data_test <- stan_model$call$sample_fn(posterior_samples_test)
post_samples_mix <- rstan::extract(fit)
print(fit)

# Extract the posterior samples for mixture (non-trunc)
alpha_samp <- post_samples_mix$alpha
mu1_samp <- post_samples_mix$mu[,1]
mu2_samp <- post_samples_mix$mu[,2]
sigma1_samp <- post_samples_mix$sigma[,1]
sigma2_samp <- post_samples_mix$sigma[,2]

#Quantiles for Table
percentiles_alpha <- quantile(alpha_samp, probs = c(0.1, 0.9))
percentiles_alpha
mean(alpha_samp)
percentiles_mu1 <- quantile(mu1_samp, probs = c(0.1, 0.9))
percentiles_mu1
mean(mu1_samp)
percentiles_mu2 <- quantile(mu2_samp, probs = c(0.1, 0.9))
percentiles_mu2
mean(mu2_samp)
percentiles_sigma1 <- quantile(sigma1_samp, probs = c(0.1, 0.9))
percentiles_sigma1
mean(sigma1_samp)
percentiles_sigma2 <- quantile(sigma2_samp, probs = c(0.1, 0.9))
percentiles_sigma2
mean(sigma2_samp)

alpha <- post_samples_mix$alpha
mu1 <- post_samples_mix$mu1
mu2 <- post_samples_mix$mu2
sigma1 <- post_samples_mix$sigma[,1]
sigma2 <- post_samples_mix$sigma[,2]

#Extract the posterior predictive samples
Y_hat_samples <- rstan::extract(fit, "Y_hat") #4000 rows of 17686 columns

#Plot hist of observed & generated data
hist(Y, main = "Observed Data")
Y_hat_test <- colMeans(post_pred_samp)
hist(Y_hat_test, main = "Generated Data")

alpha<- post_samples_mix$alpha
mu1 <- post_samples_mix$mu[,1]
mu2 <- post_samples_mix$mu[,2]
sigma1 <- post_samples_mix$sigma[,1]
sigma2 <- post_samples_mix$sigma[,2]

# Generate the sequence for x-axis
x_seq <- seq(min(Y), max(Y), length.out = 4000)

# Calculate the densities
density_mix <- rowMeans(sapply(1:length(alpha_samp), function(i) alpha_samp[i] * dnorm(x_seq, mu1_samp[i], sigma1_samp[i]) + (1 - alpha_samp[i]) * dnorm(x_seq, mu2_samp[i], sigma2_samp[i])))
density1 <- rowMeans(sapply(1:length(mu1_samp), function(i) alpha_samp[i] * dnorm(x_seq, mu1_samp[i], sigma1_samp[i])))
density2 <- rowMeans(sapply(1:length(mu2_samp), function(i) (1 - alpha_samp[i]) * dnorm(x_seq, mu2_samp[i], sigma2_samp[i])))

# Plot the histogram and the densities
hist(Y, freq = FALSE, main = "Home Run Mixture Model", xlab = "Home Run Proportion", ylab = "Density", ylim = c(0, 30))
lines(x_seq, density_mix, col = "blue")
lines(x_seq, density1, col = "red")
lines(x_seq, density2, col = "green")
legend("topright", legend = c("Mixture", "Non-Elite", "Elite"), col = c("blue", "red", "green"), lty = 1)

#Posterior Predictive Check
test_n <- length(Y) #length of data
Y_pred <- matrix(NA, nrow = test_n, ncol = 4000)

for (i in 1:4000) {
  alpha <- alpha_samp[i]
  mu1 <- mu1_samp[i]
  mu2 <- mu2_samp[i]
  sigma1 <- sigma1_samp[i]
  sigma2 <- sigma2_samp[i]
  theta <- rbinom(test_n, size = 1, prob = alpha)
  Y_pred[, i] <- rnorm(test_n, mean = (theta) * mu1 + (1-theta) * mu2, sd = (theta) * sigma1 + (1-theta) * sigma2)
}

summary_observed <- c(min = min(Y),
                      max = max(Y),
                      mean = mean(Y),
                      sd = sd(Y))

# Calculate summary statistics for each set of predicted data
summary_predicted <- apply(Y_pred, 2, function(y) {
  c(min = min(y, na.rm=TRUE),
    max = max(y, na.rm=TRUE),
    mean = mean(y, na.rm=TRUE),
    sd = sd(y, na.rm=TRUE))
})

pp_min <- summary_predicted[1,]
pp_max <- summary_predicted[2,]
pp_mean <- summary_predicted[3,]
pp_sd <- summary_predicted[4,]

ggplot(data.frame(pp_mean), aes(x = pp_mean)) +
  geom_histogram() +
  geom_vline(xintercept = mean(Y), color = "red")+
  labs(x = "Mean HR Rate")
ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram() +
  geom_vline(xintercept = min(Y), color = "red")+
  labs(x = "Minimum HR Rate")
ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram() +
  geom_vline(xintercept = max(Y), color = "red")+
  labs(x = "Max HR Rate")
ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram() +
  geom_vline(xintercept = sd(Y), color = "red")+
  labs(x = "SD HR Rate")
hist(summary_predicted[2,])
hist(summary_predicted[3,])
hist(summary_predicted[4,])


#Ranking
print(fit)
hitters <- allseason_min_ab
library(xtable)
ortiz<-hitters[hitters[,1]=="ortizda01",] 
ortiz.HRavg<-ortiz$HR/ortiz$AB
prob.elite.ortiz <- Estep(ortiz.HRavg, alpha, mu1, mu2, sigsq) 
prob.elite.ortiz <- cbind(hitters[hitters[,1]=="ortizda01",2],
                          prob.elite.ortiz)
xtable(prob.elite.ortiz , digits = c(0,0,6))

#Getting Individual probabilities for each player
Estep <- function(y, alpha, mu1, mu2, sigma1, sigma2){
  n <- length(y)  
  ind <- rep(NA,n)
  for (i in 1:n){
    prob1 <- (1-alpha)*dnorm(y[i], mean=mu1, sd=sigma1)
    prob2 <- alpha*dnorm(y[i], mean=mu2, sd=sigma2)
    ind[i] <- prob2/(prob1+prob2)
  }
  ind
}

Estep <- function(y, alpha, mu1, mu2, sigma1, sigma2){
  n <- length(y)  
  ind <- rep(NA,n)
  for (i in 1:n){
    prob1 <- (alpha)*dnorm(y[i], mean=mu1, sd=sigma1)
    prob2 <- (1-alpha)*dnorm(y[i], mean=mu2, sd=sigma2)
    ind[i] <- prob2/(prob1+prob2)
  }
  ind
}

Estep_nopower <- function(y, alpha, mu1, mu2, sigma1, sigma2){
  n <- length(y)  
  ind <- rep(NA,n)
  for (i in 1:n){
    prob1 <- (alpha)*dnorm(y[i], mean=mu1, sd=sigma1)
    prob2 <- (1-alpha)*dnorm(y[i], mean=mu2, sd=sigma2)
    ind[i] <- prob1/(prob1+prob2)
  }
  ind
}

post_pred_mean <- apply(Y_pred, 1, mean)
min(post_pred_mean)

hist(post_pred_mean)

finalindprops <- Estep_nopower(hrprop, 
                       mean(alpha_samp), 
                       mean(mu1_samp), 
                       mean(mu2_samp), 
                       mean(sigma1_samp), 
                       mean(sigma2_samp))

finalindprops <- Estep(hrprop, 
                               mean(alpha_samp), 
                               mean(mu1_samp), 
                               mean(mu2_samp), 
                               mean(sigma1_samp), 
                               mean(sigma2_samp))


hist(finalindprops)
sum(finalindprops > 0.5)
sum(finalindprops > 0.5)/length(finalindprops)
sum(finalindprops > 0.9999)
players.topHR<-allseason_min_ab[finalindprops > 0.9999,]
top_hr_hitter_2021 <- players.topHR %>% filter(yearID == 2021)

top_sluggers_2021 <- top_hr_hitter_2021$playerID

ortiz<-hitters[hitters[,1]=="ortizda01",] 
ortiz.HRavg<-ortiz$HR/ortiz$AB
prob.elite.ortiz <- Estep(ortiz.HRavg, alpha, mu1, mu2, sigsq) 
prob.elite.ortiz <- cbind(hitters[hitters[,1]=="ortizda01",2],
                          prob.elite.ortiz)


top_sluggers_2021
sort(finalindprops, decreasing =TRUE)


