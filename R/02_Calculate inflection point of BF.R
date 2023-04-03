

### Estimate polynomial regressions, find where 2nd derivative = 0, and compute separate scaling curves for each species ###

library(dplyr)
library(ggplot2)
library(car)
library(Deriv)
library(R2jags)
library(doParallel)
library(broom)
library(DHARMa)
library(AICcmodavg)

source("R/helper_functions.R")


#################
### Load data ###
#################

dat <- read_csv("Processed_data/Bite force results.csv")

cleu.bf <- dat %>%
  filter(Species == "Bull")
clim.bf <- dat %>%
  filter(Species == "Blacktip")
stib.bf <- dat %>%
  filter(Species == "Bonnethead")


ggplot(dat, aes(FL, ABF, color = Species)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ Species, scales = "free", ncol = 1)



### Fit polynomial regressions ###

dat2 <- dat %>%
  nest(.by = Species)

# Fit models
dat2 <- dat2 %>%
  mutate(poly3 = map(data, ~lm(ABF ~ poly(FL, 3, raw = T), data = .x)),
         poly4 = map(data, ~lm(ABF ~ poly(FL, 4, raw = T), data = .x))) #%>%
#   mutate(coeffs = map(mod, coefficients)) %>%
#   mutate(deriv.2 = map(coeffs, ~Deriv(f(.x), "x", nderiv = 2)))
#
# uniroot(function(x) eval(dat2$deriv.2[[2]]), c(min(stib.bf$FL),max(stib.bf$FL)))


# Which regression fits better?
AIC(dat2$poly3[[1]], dat2$poly4[[1]])  #3rd order for blacktip
AIC(dat2$poly3[[2]], dat2$poly4[[2]])  #4th order for bonnethead
AIC(dat2$poly3[[3]], dat2$poly4[[3]])  #3rd order for bull



## Viz fit of polynomials
poly.mods <- list(clim = dat2$poly3[[1]],
                  stib = dat2$poly4[[2]],
                  cleu = dat2$poly3[[3]])

fl.seq <- dat2$data %>%
  map(., ~{fl <- .x %>%
      pull(FL)

      data.frame(FL = seq(min(fl), max(fl), length.out = 100))}
      )

poly.pred <- fl.seq %>%
  bind_rows(.id = "Species") %>%
  mutate(ABF = map2(.x = poly.mods, .y = fl.seq,
                    ~predict(.x, .y)) %>%
    unlist()) %>%
  mutate(Species = case_when(Species == 1 ~ "Blacktip",
                             Species == 2 ~ "Bonnethead",
                             Species == 3 ~ "Bull"))

ggplot() +
  geom_point(data = dat, aes(FL, ABF, color = Species)) +
  geom_line(data = poly.pred, aes(FL, ABF, group = Species)) +
  theme_bw() +
  facet_wrap(~ Species, scales = "free", ncol = 1)









### Find Where 2nd Deriv = 0 ###


deriv.2 <- list(clim = coefficients(dat2$poly3[[1]]),
                stib = coefficients(dat2$poly4[[2]]),
                cleu = coefficients(dat2$poly3[[3]])) %>%
  map(., ~Deriv(f(.x), "x", nderiv = 2))

fl.range <- list(clim = c(75, 110),
                 stib = c(65, 90),
                 cleu = c(90, 180))

inflect.pt <- deriv.2[1:2] %>%
  map2(.x = ., .y = fl.range[1:2],
       ~uniroot(function(x) eval(.x), .y)$root)
# Clim = 98.8 cm FL
# Stib = 82.1

# Since root could not be found for bull sharks, I will use 120 cm FL as my threshold since this is midpoint of size range



# Viz plot w/ thresholds

ggplot() +
  geom_point(data = dat, aes(FL, ABF, color = Species)) +
  geom_line(data = poly.pred, aes(FL, ABF, group = Species)) +
  geom_vline(data = data.frame(Species = unique(dat$Species),
                               FL = c(98.8, 82.1, 120)),
             aes(xintercept = FL), color = "red") +
  theme_bw() +
  facet_wrap(~ Species, scales = "free", ncol = 1)
