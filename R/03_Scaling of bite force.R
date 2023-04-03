
### Estimate coefficients for allometric scaling of BF ###

library(tidyverse)
library(brms)


#################
### Load data ###
#################

dat <- read_csv("Processed_data/Bite force results.csv")

# Extend thresholds 5 cm either way to better estimate trend

cleu.bf <- dat %>%
  filter(Species == "Bull") %>%
  mutate(class = case_when(FL <= 120 + 5 ~ "Small",
                           FL > 120 - 5 ~ "Large"),
         log.FL = log(FL),
         log.FL.c = log.FL - mean(log.FL),
         log.ABF = log(ABF))
clim.bf <- dat %>%
  filter(Species == "Blacktip") %>%
  mutate(class = case_when(FL <= 98.8 + 5 ~ "Small",
                           FL > 98.8 - 5 ~ "Large"),
         log.FL = log(FL),
         log.FL.c = log.FL - mean(log.FL),
         log.ABF = log(ABF))
stib.bf <- dat %>%
  filter(Species == "Bonnethead") %>%
  mutate(class = case_when(FL <= 82.1 + 5 ~ "Small",
                           FL > 82.1 - 5 ~ "Large"),
         log.FL = log(FL),
         log.FL.c = log.FL - mean(log.FL),
         log.ABF = log(ABF))




####################################
### Fit allometric scaing models ###
####################################

# Define priors for scaling coeffs
prior1 <- prior(normal(5, 2.5), nlpar = "a", lb = 0) +
  prior(normal(2, 5), nlpar = "b", lb = 0) +
  prior(normal(10, 100), class = "sigma", lb = 0)



### C. leucas ###

## Small

#Prior predictive check with sample_prior="only"
fitPrior_cleu_small <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = cleu.bf %>%
    filter(class == "Small"),
  prior = prior1,
  sample_prior = "only")

#Plots the prior predictive check
pp_check(fitPrior_cleu_small, ndraws = 200) + theme_bw()



fit_cleu_small <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = cleu.bf %>%
    filter(class == "Small"),
  prior = prior1,
  iter = 4000, warmup = 2000, cores = 4, chains = 4, seed = 2023)

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
  prior = prior1,
  sample_prior = "only")

#Plots the prior predictive check
pp_check(fitPrior_cleu_large, ndraws = 200) + theme_bw()



fit_cleu_large <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = cleu.bf %>%
    filter(class == "Large"),
  prior = prior1,
  iter = 4000, warmup = 2000, cores = 4, chains = 4, seed = 2023,
  control = list(adapt_delta = 0.99))

summary(fit_cleu_large)
plot(fit_cleu_large)
pp_check(fit_cleu_large, ndraws = 200) + theme_bw()
bayes_R2(fit_cleu_large)
loo_R2(fit_cleu_large)








### C. limbatus ###

## Small

#Prior predictive check with sample_prior="only"
fitPrior_clim_small <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = clim.bf %>%
    filter(class == "Small"),
  prior = prior1,
  sample_prior = "only")

#Plots the prior predictive check
pp_check(fitPrior_clim_small, ndraws = 200) + theme_bw()



fit_clim_small <- brm(
  bf(log.ABF ~ a + log.FL.c * b, a + b ~ 1, nl = TRUE),
  family = gaussian,
  data = clim.bf %>%
    filter(class == "Small"),
  prior = prior1,
  iter = 4000, warmup = 2000, cores = 4, chains = 4, seed = 2023)

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
  iter = 4000, warmup = 1000, cores = 4, chains = 4, seed = 2023,
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
  iter = 4000, warmup = 1000, cores = 4, chains = 4, seed = 2023)

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
  iter = 4000, warmup = 1000, cores = 4, chains = 4, seed = 2023,
  control = list(adapt_delta = 0.9))

summary(fit_stib_large)
plot(fit_stib_large)
pp_check(fit_stib_large, ndraws = 200) + theme_bw()
bayes_R2(fit_stib_large)
loo_R2(fit_stib_large)
