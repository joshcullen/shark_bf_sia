
### Estimate coefficients for allometric scaling of BF ###

library(tidyverse)
library(ggdist)
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
                   seed = 2023,
                   control = list(adapt_delta = 0.9))


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
  filter(FL <= 82.1 + 1)

stib.l <- stib %>%
  filter(FL > 82.1 - 1)



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




### Combine predictions from small and large bonnethead sharks ###

ggplot() +
  ggdist::stat_lineribbon(data = y_hat_stib_s, aes(x = FL, y = exp(prediction)),
                          .width = c(.95, .5), color = "#08519C", alpha = 0.7) +
  ggdist::stat_lineribbon(data = y_hat_stib_l, aes(x = FL, y = exp(prediction)),
                          .width = c(.95, .5), color = "#08519C", alpha = 0.7) +
  geom_point(data = stib, aes(x = FL, y = ABF), size = 2) +
  scale_fill_brewer() +
  labs(y = "ABF (N)") +
  theme_bw()











