
### Estimate coefficients for allometric scaling of BF ###

library(tidyverse)
library(ggdist)
library(rstan)
library(rstantools)
library(bayesplot)
library(MetBrewer)
library(ggtext)
library(patchwork)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



#################
### Load data ###
#################

dat <- read_csv("Processed_data/Bite force results.csv")


inflect.pts <- c(blacktip = 103.5, bonnethead = 79.3)



# Separate by species and transform FL and ABF for efficient MCMC sampling

cleu <- dat %>%
  filter(Species == "Bull") %>%
  mutate(
    log.FL = log(FL),
    log.FL.c = log.FL - mean(log.FL))

clim <- dat %>%
  filter(Species == "Blacktip") %>%
  mutate(
    log.FL = log(FL),
    log.FL.c = log.FL - mean(log.FL))

stib <- dat %>%
  filter(Species == "Bonnethead") %>%
  mutate(
    log.FL = log(FL),
    log.FL.c = log.FL - mean(log.FL))


# Define color palette for viz
pal <- met.brewer(palette_name = "Egypt", n = 3, return_hex = TRUE)





#########################
### Bull shark models ###
#########################

data.cleu <- list(n = nrow(cleu),
                  ABF = cleu$ABF,
                  log_FL_c = cleu$log.FL.c)



### Fit bull shark model ###

fit.cleu <- stan(file = "R/OLS allometric scaling.stan",
                 data = data.cleu,
                 iter = 4000,
                 thin = 1,
                 chains = 4,
                 seed = 2024)


## Check model diagnostics
print(fit.cleu,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

# Traceplot
traceplot(fit.cleu, pars = c("log_a", "b", "sigma")) +
  scale_color_viridis_d(option = "mako")

# Point estimate and CI
plot(fit.cleu, pars = c("log_a", "b", "sigma"))

# Density plot
plot(fit.cleu, pars = c("log_a", "b", "sigma"),
     plotfun = "dens", separate_chains = TRUE, alpha = 0.3) +
  scale_fill_viridis_d(option = "mako")



## Posterior prediction
y_hat <- extract(fit.cleu, par = 'y_pp')$y_pp
ppc_dens_overlay(y = cleu$ABF,
                 yrep = y_hat[1:50,]) +
  theme_bw()

## Calculate R^2
R2 <- bayes_R2(y_hat, cleu$ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)  #0.88


y_hat_cleu <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = cleu$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


cleu_plot <- ggplot() +
  ggdist::stat_lineribbon(data = y_hat_cleu, aes(x = FL, y = prediction),
                  .width = c(.95, .5), color = pal[1], alpha = 0.25,
                  fill = pal[1], linewidth = 0.5) +
  geom_point(data = cleu, aes(x = FL, y = ABF), size = 2, color = pal[1], alpha = 0.7) +
  labs(x = 'FL (cm)', y = 'ABF (N)', title = "**<span style = 'color: #dd5129'>Bull</span>**") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        plot.margin=unit(rep(1, 4),'lines'),
        plot.title = element_markdown(size = 12))








#############################
### Blacktip shark models ###
#############################

clim.s <- clim %>%
  filter(FL <= inflect.pts[1] + 3)

clim.l <- clim %>%
  filter(FL >= inflect.pts[1] - 3)



data.clim.s <- list(n = nrow(clim.s),
                    ABF = clim.s$ABF,
                    log_FL_c = clim.s$log.FL.c)

data.clim.l <- list(n = nrow(clim.l),
                    ABF = clim.l$ABF,
                    log_FL_c = clim.l$log.FL.c)



### Fit small blacktip shark model ###

fit.clim.s <- stan(file = "R/OLS allometric scaling.stan",
                   data = data.clim.s,
                   iter = 4000,
                   thin = 1,
                   chains = 4,
                   seed = 2024)


## Check model diagnostics
print(fit.clim.s,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

# Traceplot
traceplot(fit.clim.s, pars = c("log_a", "b", "sigma")) +
  scale_color_viridis_d(option = "mako")

# Point estimate and CI
plot(fit.clim.s, pars = c("log_a", "b", "sigma"))

# Density plot
plot(fit.clim.s, pars = c("log_a", "b", "sigma"),
     plotfun = "dens", separate_chains = TRUE, alpha = 0.3) +
  scale_fill_viridis_d(option = "mako")


## Posterior prediction
y_hat <- extract(fit.clim.s, par = 'y_pp')$y_pp
ppc_dens_overlay(y = clim.s$ABF,
                 yrep = y_hat[1:50,]) +
  theme_bw()

## Calculate R^2
R2 <- bayes_R2(y_hat, clim.s$ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)  #0.85


y_hat_clim_s <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = clim.s$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_clim_s, aes(x = FL, y = prediction),
                          .width = c(.95, .5), color = pal[2], alpha = 0.25,
                          fill = pal[2], linewidth = 2) +
  geom_point(data = clim.s, aes(x = FL, y = ABF), size = 2) +
  labs(x = 'FL (cm)', y = 'ABF (N)', title = "Blacktip - small") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())




### Fit large blacktip shark model ###

fit.clim.l <- stan(file = "R/OLS allometric scaling.stan",
                   data = data.clim.l,
                   iter = 4000,
                   thin = 1,
                   chains = 4,
                   seed = 2024)


## Check model diagnostics
print(fit.clim.l,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

# Traceplot
traceplot(fit.clim.l, pars = c("log_a", "b", "sigma")) +
  scale_color_viridis_d(option = "mako")

# Point estimate and CI
plot(fit.clim.l, pars = c("log_a", "b", "sigma"))

# Density plot
plot(fit.clim.l, pars = c("log_a", "b", "sigma"),
     plotfun = "dens", separate_chains = TRUE, alpha = 0.3) +
  scale_fill_viridis_d(option = "mako")


## Posterior prediction
y_hat <- extract(fit.clim.l, par = 'y_pp')$y_pp
ppc_dens_overlay(y = clim.l$ABF,
                 yrep = y_hat[1:50,]) +
  theme_bw()

## Calculate R^2
R2 <- bayes_R2(y_hat, clim.l$ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)  #0.34


y_hat_clim_l <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = clim.l$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_clim_l, aes(x = FL, y = prediction),
                          .width = c(.95, .5), color = pal[2], alpha = 0.25,
                          fill = pal[2], linewidth = 2) +
  geom_point(data = clim.l, aes(x = FL, y = ABF), size = 2) +
  labs(x = 'FL (cm)', y = 'ABF (N)', title = "Blacktip - large") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())




### Combine predictions from small and large blacktip sharks ###

clim_plot <- ggplot() +
  ggdist::stat_lineribbon(data = y_hat_clim_s, aes(x = FL, y = prediction),
                          .width = c(.95, .5), color = pal[2], alpha = 0.25,
                          fill = pal[2], linewidth = 0.5) +
  ggdist::stat_lineribbon(data = y_hat_clim_l, aes(x = FL, y = prediction),
                          .width = c(.95, .5), color = pal[2], alpha = 0.25,
                          fill = pal[2], linewidth = 0.5) +
  geom_point(data = clim, aes(x = FL, y = ABF), size = 2, color = pal[2], alpha = 0.7) +
  geom_vline(xintercept = inflect.pts[1], linetype = "dashed") +
  labs(x = 'FL (cm)', y = 'ABF (N)', title = "**<span style = 'color: #0f7ba2'>Blacktip</span>**") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        plot.title = element_markdown(size = 12))









###############################
### Bonnethead shark models ###
###############################

stib.s <- stib %>%
  filter(FL <= 79.3 + 3)

stib.l <- stib %>%
  filter(FL >= 79.3 - 3)



data.stib.s <- list(n = nrow(stib.s),
                    ABF = stib.s$ABF,
                    log_FL_c = stib.s$log.FL.c)

data.stib.l <- list(n = nrow(stib.l),
                    ABF = stib.l$ABF,
                    log_FL_c = stib.l$log.FL.c)



### Fit small bonnethead shark model ###

fit.stib.s <- stan(file = "R/OLS allometric scaling.stan",
                   data = data.stib.s,
                   iter = 4000,
                   thin = 1,
                   chains = 4,
                   seed = 2024)


## Check model diagnostics
print(fit.stib.s,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

# Traceplot
traceplot(fit.stib.s, pars = c("log_a", "b", "sigma")) +
  scale_color_viridis_d(option = "mako")

# Point estimate and CI
plot(fit.stib.s, pars = c("log_a", "b", "sigma"))

# Density plot
plot(fit.stib.s, pars = c("log_a", "b", "sigma"),
     plotfun = "dens", separate_chains = TRUE, alpha = 0.3) +
  scale_fill_viridis_d(option = "mako")


## Posterior prediction
y_hat <- extract(fit.stib.s, par = 'y_pp')$y_pp
ppc_dens_overlay(y = stib.s$ABF,
                 yrep = y_hat[1:50,]) +
  theme_bw()

## Calculate R^2
R2 <- bayes_R2(y_hat, stib.s$ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)  #0.48


y_hat_stib_s <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = stib.s$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_stib_s, aes(x = FL, y = prediction),
                          .width = c(.95, .5), color = pal[3], alpha = 0.25,
                          fill = pal[3], linewidth = 2) +
  geom_point(data = stib.s, aes(x = FL, y = ABF), size = 2) +
  labs(x = 'FL (cm)', y = 'ABF (N)', title = "Bonnethead - small") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())




### Fit large bonnethead shark model ###

fit.stib.l <- stan(file = "R/OLS allometric scaling.stan",
                   data = data.stib.l,
                   iter = 4000,
                   thin = 1,
                   chains = 4,
                   seed = 2024)


## Check model diagnostics
print(fit.stib.l,
      digits_summary = 3,
      pars = c("log_a", "b", "sigma"),
      probs = c(0.025, 0.5, 0.975))

# Traceplot
traceplot(fit.stib.l, pars = c("log_a", "b", "sigma")) +
  scale_color_viridis_d(option = "mako")

# Point estimate and CI
plot(fit.stib.l, pars = c("log_a", "b", "sigma"))

# Density plot
plot(fit.stib.l, pars = c("log_a", "b", "sigma"),
     plotfun = "dens", separate_chains = TRUE, alpha = 0.3) +
  scale_fill_viridis_d(option = "mako")


## Posterior prediction
y_hat <- extract(fit.stib.l, par = 'y_pp')$y_pp
ppc_dens_overlay(y = stib.l$ABF,
                 yrep = y_hat[1:50,]) +
  theme_bw()

## Calculate R^2
R2 <- bayes_R2(y_hat, stib.l$ABF)
quantile(R2, c(0.025, 0.5, 0.95))
mean(R2)  #0.32


y_hat_stib_l <- y_hat %>%
  t() %>%
  data.frame() %>%
  mutate(FL = stib.l$FL) %>%
  pivot_longer(cols = -FL, names_to = "iter", values_to = "prediction") %>%
  group_by(FL, iter)


ggplot() +
  ggdist::stat_lineribbon(data = y_hat_stib_l, aes(x = FL, y = prediction),
                          .width = c(.95, .5), color = pal[3], alpha = 0.25,
                          fill = pal[3], linewidth = 2) +
  geom_point(data = stib.l, aes(x = FL, y = ABF), size = 2) +
  labs(x = 'FL (cm)', y = 'ABF (N)', title = "Bonnethead - large") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())




### Combine predictions from small and large bonnethead sharks ###

stib_plot <- ggplot() +
  ggdist::stat_lineribbon(data = y_hat_stib_s, aes(x = FL, y = prediction),
                          .width = c(.95, .5), color = pal[3], alpha = 0.25,
                          fill = pal[3], linewidth = 0.5) +
  ggdist::stat_lineribbon(data = y_hat_stib_l, aes(x = FL, y = prediction),
                          .width = c(.95, .5), color = pal[3], alpha = 0.25,
                          fill = pal[3], linewidth = 0.5) +
  geom_point(data = stib, aes(x = FL, y = ABF), size = 2, color = pal[3], alpha = 0.7) +
  geom_vline(xintercept = inflect.pts[2], linetype = "dashed") +
  labs(x = 'FL (cm)', y = 'ABF (N)', title = "**<span style = 'color: #43b284'>Bonnethead</span>**") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        plot.title = element_markdown(size = 12))





### Combine plots from all species ###

cleu_plot + clim_plot + stib_plot +
  plot_layout(nrow = 3) &
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggsave("Figures/ABF scaling estimates.png", width = 5, height = 7, units = "in", dpi = 600)




