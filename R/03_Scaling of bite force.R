
### Estimate coefficients for allometric scaling of BF ###

library(tidyverse)
# library(brms)
library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



#################
### Load data ###
#################

dat <- read_csv("Processed_data/Bite force results.csv")

# Extend thresholds 10 cm either way to better estimate trend

cleu <- dat %>%
  filter(Species == "Bull") %>%
  mutate(
    # class = case_when(FL <= 120 + 5 ~ "Small",
    #                   FL > 120 - 5 ~ "Large"),
    log.FL = log(FL),
    log.FL.c = log.FL - mean(log.FL),
    log.ABF = log(ABF))

clim <- dat %>%
  filter(Species == "Blacktip") %>%
  mutate(
    # class = case_when(FL <= 98.8 + 10 ~ "Small",
    #                   FL > 98.8 - 10 ~ "Large"),
    log.FL = log(FL),
    log.FL.c = log.FL - mean(log.FL),
    log.ABF = log(ABF))

stib <- dat %>%
  filter(Species == "Bonnethead") %>%
  mutate(
    # class = case_when(FL <= 82.1 + 5 ~ "Small",
    #                   FL > 82.1 - 5 ~ "Large"),
    log.FL = log(FL),
    log.FL.c = log.FL - mean(log.FL),
    log.ABF = log(ABF))




#########################
### Bull shark models ###
#########################

cleu.s <- cleu %>%
  filter(FL <= 120 + 5)

cleu.l <- cleu %>%
  filter(FL >= 120 - 5)



data.cleu.s <- list(n = nrow(cleu.s),
            log_ABF = cleu.s$log.ABF,
            log_FL_c = cleu.s$log.FL.c)

data.cleu.l <- list(n = nrow(cleu.l),
                    log_ABF = cleu.l$log.ABF,
                    log_FL_c = cleu.l$log.FL.c)



### Fit small bull shark model ###

fit.cleu.s <- stan(file = "R/Stan MA regression.stan",
                  data = data.cleu.s,
                  iter = 3000,
                  warmup = 1000,
                  thin = 1,
                  chains = 4,
                  seed = 2023)


## Check model diagnostics
print(fit.cleu.s,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

mcmc_intervals(fit.cleu.s,
           pars = c("log_a", "b", "sigma"),
           prob = 0.95, prob_outer = 0.99) +
  theme_bw()

mcmc_dens_overlay(fit.cleu.s, pars = c("log_a", "b", "sigma")) +
  theme_bw()

mcmc_trace(fit.cleu.s, pars = c("log_a", "b", "sigma")) +
  theme_bw()


## Posterior prediction
y_hat <- extract(fit.cleu.s, par = 'y_pp')$y_pp

ppc_dens_overlay(y = exp(cleu.s$log.ABF),
                 yrep = exp(y_hat[1:50,])) +
  theme_bw()

## Calculate R^2
R2 <- rstantools::bayes_R2(y_hat, cleu.s$log.ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)


y_hat_cleu_s <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = cleu.s$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_cleu_s, aes(x = FL, y = exp(prediction)),
                  .width = c(.99, .95, .5), color = "#08519C") +
  geom_point(data = cleu.s, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  theme_bw()




### Fit large bull shark model ###

fit.cleu.l <- stan(file = "R/Stan MA regression.stan",
                   data = data.cleu.l,
                   iter = 3000,
                   warmup = 1000,
                   thin = 1,
                   chains = 4,
                   seed = 2023)


## Check model diagnostics
print(fit.cleu.l,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

mcmc_intervals(fit.cleu.l,
               pars = c("log_a", "b", "sigma"),
               prob = 0.95, prob_outer = 0.99) +
  theme_bw()

mcmc_dens_overlay(fit.cleu.l, pars = c("log_a", "b", "sigma")) +
  theme_bw()

mcmc_trace(fit.cleu.l, pars = c("log_a", "b", "sigma")) +
  theme_bw()


## Posterior prediction
y_hat <- extract(fit.cleu.l, par = 'y_pp')$y_pp

ppc_dens_overlay(y = exp(cleu.l$log.ABF),
                 yrep = exp(y_hat[1:50,])) +
  theme_bw()

## Calculate R^2
R2 <- rstantools::bayes_R2(y_hat, cleu.l$log.ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)


y_hat_cleu_l <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = cleu.l$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_cleu_l, aes(x = FL, y = exp(prediction)),
                          .width = c(.99, .95, .5), color = "#08519C") +
  geom_point(data = cleu.l, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  theme_bw()




### Combine predictions from small and large bull sharks ###

ggplot() +
  ggdist::stat_lineribbon(data = y_hat_cleu_s, aes(x = FL, y = exp(prediction)),
                          .width = c(.95, .5), color = "#08519C", alpha = 0.7) +
  ggdist::stat_lineribbon(data = y_hat_cleu_l, aes(x = FL, y = exp(prediction)),
                          .width = c(.95, .5), color = "#08519C", alpha = 0.7) +
  geom_point(data = cleu, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  labs(y = "ABF (N)") +
  theme_bw()









#############################
### Blacktip shark models ###
#############################

clim.s <- clim %>%
  filter(FL <= 98.8 + 5)

clim.l <- clim %>%
  filter(FL >= 98.8 - 5)



data.clim.s <- list(n = nrow(clim.s),
                    log_ABF = clim.s$log.ABF,
                    log_FL_c = clim.s$log.FL.c)

data.clim.l <- list(n = nrow(clim.l),
                    log_ABF = clim.l$log.ABF,
                    log_FL_c = clim.l$log.FL.c)



### Fit small blacktip shark model ###

fit.clim.s <- stan(file = "R/Stan MA regression.stan",
                   data = data.clim.s,
                   iter = 3000,
                   warmup = 1000,
                   thin = 1,
                   chains = 4,
                   seed = 2023)


## Check model diagnostics
print(fit.clim.s,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

mcmc_intervals(fit.clim.s,
               pars = c("log_a", "b", "sigma"),
               prob = 0.95, prob_outer = 0.99) +
  theme_bw()

mcmc_dens_overlay(fit.clim.s, pars = c("log_a", "b", "sigma")) +
  theme_bw()

mcmc_trace(fit.clim.s, pars = c("log_a", "b", "sigma")) +
  theme_bw()


## Posterior prediction
y_hat <- extract(fit.clim.s, par = 'y_pp')$y_pp

ppc_dens_overlay(y = exp(clim.s$log.ABF),
                 yrep = exp(y_hat[1:50,])) +
  theme_bw()

## Calculate R^2
R2 <- rstantools::bayes_R2(y_hat, clim.s$log.ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)


y_hat_clim_s <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = clim.s$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_clim_s, aes(x = FL, y = exp(prediction)),
                          .width = c(.99, .95, .5), color = "#08519C") +
  geom_point(data = clim.s, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  theme_bw()




### Fit large blacktip shark model ###

fit.clim.l <- stan(file = "R/Stan MA regression.stan",
                   data = data.clim.l,
                   iter = 3000,
                   warmup = 1000,
                   thin = 1,
                   chains = 4,
                   seed = 2023)


## Check model diagnostics
print(fit.clim.l,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

mcmc_intervals(fit.clim.l,
               pars = c("log_a", "b", "sigma"),
               prob = 0.95, prob_outer = 0.99) +
  theme_bw()

mcmc_dens_overlay(fit.clim.l, pars = c("log_a", "b", "sigma")) +
  theme_bw()

mcmc_trace(fit.clim.l, pars = c("log_a", "b", "sigma")) +
  theme_bw()


## Posterior prediction
y_hat <- extract(fit.clim.l, par = 'y_pp')$y_pp

ppc_dens_overlay(y = exp(clim.l$log.ABF),
                 yrep = exp(y_hat[1:50,])) +
  theme_bw()

## Calculate R^2
R2 <- rstantools::bayes_R2(y_hat, clim.l$log.ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)


y_hat_clim_l <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = clim.l$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_clim_l, aes(x = FL, y = exp(prediction)),
                          .width = c(.99, .95, .5), color = "#08519C") +
  geom_point(data = clim.l, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  theme_bw()




### Combine predictions from small and large blacktip sharks ###

ggplot() +
  ggdist::stat_lineribbon(data = y_hat_clim_s, aes(x = FL, y = exp(prediction)),
                          .width = c(.95, .5), color = "#08519C", alpha = 0.7) +
  ggdist::stat_lineribbon(data = y_hat_clim_l, aes(x = FL, y = exp(prediction)),
                          .width = c(.95, .5), color = "#08519C", alpha = 0.7) +
  geom_point(data = clim, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  labs(y = "ABF (N)") +
  theme_bw()









###############################
### Bonnethead shark models ###
###############################

stib.s <- stib %>%
  filter(FL <= 82.1 + 0)

stib.l <- stib %>%
  filter(FL > 82.1 - 0)



data.stib.s <- list(n = nrow(stib.s),
                    log_ABF = stib.s$log.ABF,
                    log_FL_c = stib.s$log.FL.c)

data.stib.l <- list(n = nrow(stib.l),
                    log_ABF = stib.l$log.ABF,
                    log_FL_c = stib.l$log.FL.c)



### Fit small bonnethead shark model ###

fit.stib.s <- stan(file = "R/Stan MA regression.stan",
                   data = data.stib.s,
                   iter = 3000,
                   warmup = 1000,
                   thin = 1,
                   chains = 4,
                   seed = 2023)


## Check model diagnostics
print(fit.stib.s,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

mcmc_intervals(fit.stib.s,
               pars = c("log_a", "b", "sigma"),
               prob = 0.95, prob_outer = 0.99) +
  theme_bw()

mcmc_dens_overlay(fit.stib.s, pars = c("log_a", "b", "sigma")) +
  theme_bw()

mcmc_trace(fit.stib.s, pars = c("log_a", "b", "sigma")) +
  theme_bw()


## Posterior prediction
y_hat <- extract(fit.stib.s, par = 'y_pp')$y_pp

ppc_dens_overlay(y = exp(stib.s$log.ABF),
                 yrep = exp(y_hat[1:50,])) +
  theme_bw()

## Calculate R^2
R2 <- rstantools::bayes_R2(y_hat, stib.s$log.ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)


y_hat_stib_s <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = stib.s$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_stib_s, aes(x = FL, y = exp(prediction)),
                          .width = c(.99, .95, .5), color = "#08519C") +
  geom_point(data = stib.s, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  theme_bw()




### Fit large bonnethead shark model ###

fit.stib.l <- stan(file = "R/Stan MA regression.stan",
                   data = data.stib.l,
                   iter = 3000,
                   warmup = 1000,
                   thin = 1,
                   chains = 4,
                   seed = 2023)


## Check model diagnostics
print(fit.stib.l,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

mcmc_intervals(fit.stib.l,
               pars = c("log_a", "b", "sigma"),
               prob = 0.95, prob_outer = 0.99) +
  theme_bw()

mcmc_dens_overlay(fit.stib.l, pars = c("log_a", "b", "sigma")) +
  theme_bw()

mcmc_trace(fit.stib.l, pars = c("log_a", "b", "sigma")) +
  theme_bw()


## Posterior prediction
y_hat <- extract(fit.stib.l, par = 'y_pp')$y_pp

ppc_dens_overlay(y = exp(stib.l$log.ABF),
                 yrep = exp(y_hat[1:50,])) +
  theme_bw()

## Calculate R^2
R2 <- rstantools::bayes_R2(y_hat, stib.l$log.ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)


y_hat_stib_l <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = stib.l$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_stib_l, aes(x = FL, y = exp(prediction)),
                          .width = c(.99, .95, .5), color = "#08519C") +
  geom_point(data = stib.l, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  theme_bw()




### Combine predictions from small and large blacktip sharks ###

ggplot() +
  ggdist::stat_lineribbon(data = y_hat_stib_s, aes(x = FL, y = exp(prediction)),
                          .width = c(.95, .5), color = "#08519C", alpha = 0.7) +
  ggdist::stat_lineribbon(data = y_hat_stib_l, aes(x = FL, y = exp(prediction)),
                          .width = c(.95, .5), color = "#08519C", alpha = 0.7) +
  geom_point(data = stib, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  labs(y = "ABF (N)") +
  theme_bw()

















####################################
### Fit allometric scaling models ###
####################################

# Define priors for scaling coeffs
prior1 <- prior(normal(5, 2.5), nlpar = "a", lb = 0) +
  prior(normal(2, 1), nlpar = "b", lb = 0) +
  prior(normal(2, 2), class = "sigma", lb = 0)



### C. leucas ###

## Small

#Prior predictive check with sample_prior="only"
fitPrior_cleu_small <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = cleu.bf %>%
    filter(class == "Small"),
  prior = prior(normal(5, 2.5), nlpar = "a", lb = 0) +
    prior(normal(2, 5), nlpar = "b", lb = 0) +
    prior(normal(10, 10), class = "sigma", lb = 0),
  sample_prior = "only")

#Plots the prior predictive check
pp_check(fitPrior_cleu_small, ndraws = 200) + theme_bw()



fit_cleu_small <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = cleu.bf %>%
    filter(class == "Small"),
  prior = prior(normal(5, 2.5), nlpar = "a", lb = 0) +
    prior(normal(2, 5), nlpar = "b", lb = 0) +
    prior(normal(10, 10), class = "sigma", lb = 0),
  iter = 7000, warmup = 2000, cores = 4, chains = 4, seed = 2023)

summary(fit_cleu_small)
plot(fit_cleu_small)
pp_check(fit_cleu_small, ndraws = 200) + theme_bw()
bayes_R2(fit_cleu_small)
loo_R2(fit_cleu_small)



## Large

#Prior predictive check with sample_prior="only"
fitPrior_cleu_large <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = cleu.bf %>%
    filter(class == "Large"),
  prior = prior(normal(5, 2.5), nlpar = "a", lb = 0) +
    prior(normal(2, 5), nlpar = "b", lb = 0) +
    prior(normal(10, 10), class = "sigma", lb = 0),
  sample_prior = "only")

#Plots the prior predictive check
pp_check(fitPrior_cleu_large, ndraws = 200) + theme_bw()



fit_cleu_large <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = cleu.bf %>%
    filter(class == "Large"),
  prior = prior1,
  iter = 7000, warmup = 2000, cores = 4, chains = 4, seed = 2023,
  control = list(adapt_delta = 0.99))

summary(fit_cleu_large)
plot(fit_cleu_large)
pp_check(fit_cleu_large, ndraws = 200) + theme_bw()
bayes_R2(fit_cleu_large)
loo_R2(fit_cleu_large)








### C. limbatus ###

## Small

clim.bf.s <- clim.bf %>%
  filter(FL <= 98.8 + 10)

#Prior predictive check with sample_prior="only"
fitPrior_clim_small <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = clim.bf.s,
  prior = prior(normal(5, 100), nlpar = "a", lb = 0) +
    prior(normal(2, 100), nlpar = "b", lb = 0) +
    prior(gamma(1, 1), class = "sigma"),
  sample_prior = "only")

#Plots the prior predictive check
pp_check(fitPrior_clim_small, ndraws = 200) + theme_bw()



fit_clim_small <- brm(
  # bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  log.ABF | mi(log.ABF.sd) ~ me(log.FL.c, log.FL.sd),
  family = gaussian,
  data = clim.bf.s,
  prior = c(prior(normal(5, 100), class = "Intercept", lb = 0),
    prior(normal(2, 100), class = "b", lb = 0),
    prior(normal(0, 1), class = "meanme"),
    prior(gamma(1, 1), class = "sigma")),
  backend = 'cmdstanr',
  iter = 7000, warmup = 2000, cores = 4, chains = 4, seed = 2023)

summary(fit_clim_small)
plot(fit_clim_small)
pp_check(fit_clim_small, ndraws = 200) + theme_bw()
bayes_R2(fit_clim_small)
loo_R2(fit_clim_small)



## Large

#Prior predictive check with sample_prior="only"
fitPrior_clim_large <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = clim.bf %>%
    filter(class == "Large"),
  prior = prior1,
  sample_prior = "only")

#Plots the prior predictive check
pp_check(fitPrior_clim_large, ndraws = 200) + theme_bw()



fit_clim_large <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = clim.bf %>%
    filter(class == "Large"),
  prior = prior1,
  iter = 7000, warmup = 1000, cores = 4, chains = 4, seed = 2023,
  control = list(adapt_delta = 0.9))

summary(fit_clim_large)
plot(fit_clim_large)
pp_check(fit_clim_large, ndraws = 200) + theme_bw()
bayes_R2(fit_clim_large)
loo_R2(fit_clim_large)






### S. tiburo ###

## Small

#Prior predictive check with sample_prior="only"
fitPrior_stib_small <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = stib.bf %>%
    filter(class == "Small"),
  prior = prior1,
  sample_prior = "only")

#Plots the prior predictive check
pp_check(fitPrior_stib_small, ndraws = 200) + theme_bw()



fit_stib_small <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = stib.bf %>%
    filter(class == "Small"),
  prior = prior1,
  iter = 7000, warmup = 1000, cores = 4, chains = 4, seed = 2023)

summary(fit_stib_small)
plot(fit_stib_small)
pp_check(fit_stib_small, ndraws = 200) + theme_bw()
bayes_R2(fit_stib_small)
loo_R2(fit_stib_small)



## Large

#Prior predictive check with sample_prior="only"
fitPrior_stib_large <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = stib.bf %>%
    filter(class == "Large"),
  prior = prior1,
  sample_prior = "only")

#Plots the prior predictive check
pp_check(fitPrior_stib_large, ndraws = 200) + theme_bw()



fit_stib_large <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = stib.bf %>%
    filter(class == "Large"),
  prior = prior1,
  iter = 7000, warmup = 1000, cores = 4, chains = 4, seed = 2023,
  control = list(adapt_delta = 0.99))

summary(fit_stib_large)
plot(fit_stib_large)
pp_check(fit_stib_large, ndraws = 200) + theme_bw()
bayes_R2(fit_stib_large)
loo_R2(fit_stib_large)



####################
### JAGS Version ###
####################

library(R2jags)
library(tidybayes)
library(broom.mixed)

##Data
ABF<-clim.bf %>%
  filter(class == "Small") %>%
  pull(ABF)
log.ABF <- log(ABF)
FL<- clim.bf %>%
  filter(class == "Small") %>%
  pull(FL)
log.FL<- log(FL)
log.FL.c<- log.FL - mean(log.FL)
n<- clim.bf %>%
  filter(class == "Small") %>%
  nrow()
# clim.small2$log.FL.c<- (log(clim.small2$FL) - mean(log(clim.small2$FL)))

clim.small.data<- list("log.ABF","log.FL.c","n")

##Model

clim.small.model.c2<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log.a + b*log.FL.c[k]
    # bf[k]<- exp(log.bf[k])
    log.ABF[k] ~ dnorm(log.bf[k], tau)
  }
  #Priors
  # log.a ~ dunif(0,8) #changed from 50 to 500 to raise ceiling
  # b ~ dunif(0.01,15)
  # sigma ~ dunif(0,45*2)
  log.a ~ dnorm(5,0.001) #changed from 50 to 500 to raise ceiling
  b ~ dnorm(2,0.001)
  sigma ~ dgamma(1,1)
  tau <- 1/sigma/sigma
  # Posterior predictive simulations
  for (i in 1:n) {
    log.bf.post[i] <- log.a + b*log.FL.c[i]
    log.abfPred[i] ~ dnorm(log.bf.post[i], tau)
    abfPred[i] <- exp(log.abfPred[i])
  }
}


##Params
clim.small.params2<- c("log.a","b","sigma","abfPred")


##Inits

clim.small.inits2<- function() {
  #for 'a', 50 changed to 500 to raise ceiling
  list("log.a" = rnorm(1,5,10), "b" = rnorm(1,2,10), "sigma" = rgamma(1,1,1))
}

##Run model

clim.small.fit.c4<-jags(data=clim.small.data,
                        inits=clim.small.inits2,
                        parameters.to.save=clim.small.params2,
                        model.file=clim.small.model.c2,
                        n.iter=10000*10+100, n.burnin=100, n.thin=10, n.chains=3)
clim.small.mcmc.c4<- as.mcmc(clim.small.fit.c4)

clim.small.fit.c4

plot(clim.small.mcmc.c4[,c("log.a","b","sigma")]) #mixing looks good
acfplot(clim.small.mcmc.c4[,c("log.a","b","sigma")]) #looks good for 'a' and 'b'

#Run diagnostic tests
gelman.diag(clim.small.mcmc.c4[,c("log.a","b","sigma")]); gelman.plot(clim.small.mcmc.c4[,c("log.a","b","sigma")]) #looks very good


simulations = clim.small.fit.c4$BUGSoutput$sims.list$abfPred

sim.df <- data.frame(simulations) %>%
  mutate(iter = 1:n()) %>%
  slice_sample(n = 200) %>%
  pivot_longer(cols = -iter, names_to = 'ID', values_to = 'abfPred')

ggplot() +
  geom_density(data = sim.df, aes(abfPred, group = iter), color = alpha("lightblue", 0.3)) +
  geom_density(data = clim.bf %>%
                 filter(class == 'Small'), aes(ABF)) +
  theme_bw()





#create new power fxn

a<- exp(as.numeric(clim.small.fit.c4$BUGSoutput$mean$log.a))
modifier.for.a <- exp(-b*mean(log(FL)))
b<- as.numeric(clim.small.fit.c4$BUGSoutput$mean$b)

clim.small.fxn.alt <- function(x) a*modifier.for.a*x^b

ggplot(data = data.frame(FL = FL,
                         ABF = ABF), aes(x = FL, y = ABF)) +
  geom_point(size = 3) +
  stat_function(fun = clim.small.fxn.alt, size = 1, color = "blue") +
  theme_bw()


#Bayesian R2 (via Gelman et al 2017)

mcmc.clim.small.fit <- clim.small.fit.c4$BUGSoutput$sims.matrix
Xmat.clim.small.fit = model.matrix(~log.FL.c, clim.bf %>%
                                     filter(class == 'Small')
                                   )
coefs.clim.small.fit = mcmc.clim.small.fit[, c("log.a", "b")]
fit_clim.small.fit = coefs.clim.small.fit %*% t(Xmat.clim.small.fit)
resid_clim.small.fit = sweep(fit_clim.small.fit, 2, log(ABF), "-")
var_f_clim.small = apply(fit_clim.small.fit, 1, var)
var_e_clim.small = apply(resid_clim.small.fit, 1, var)
R2_clim.small = var_f_clim.small/(var_f_clim.small + var_e_clim.small)
tidyMCMC(as.mcmc(R2_clim.small), conf.int = TRUE, conf.method = "HPDinterval")






clim.bf.s %>%
  modelr::data_grid(log.FL.c = modelr::seq_range(log.FL.c, n = 101)) %>%
  add_predicted_draws(fit_clim_small) %>%
  ggplot(aes(x = log.FL.c, y = exp(log.ABF))) +
  stat_lineribbon(aes(y = exp(.prediction)), .width = c(.99, .95, .9, .5), color = "#08519C") +
  geom_point(data = clim.bf.s, size = 2) +
  scale_fill_brewer() +
  theme_bw()


foo <- clim.small.fit.c4$BUGSoutput$sims.list$abfPred %>%
  t() %>%
  data.frame() %>%
  mutate(FL = FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  stat_lineribbon(data = foo, aes(x = FL, y = prediction),
                  .width = c(.99, .95, .5), color = "#08519C") +
  geom_point(data = clim.bf %>%
               filter(class == 'Small'), aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  theme_bw()
