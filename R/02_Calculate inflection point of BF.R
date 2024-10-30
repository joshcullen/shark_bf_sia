

### Estimate inflection point using derivatives of GAMs ###

library(tidyverse)
library(mgcv)
library(gratia)
library(MetBrewer)



#################
### Load data ###
#################

dat <- read_csv("Processed_data/Bite force results.csv") |>
  mutate(Species = factor(Species, levels = c("Bull","Blacktip","Bonnethead")))


ggplot(dat, aes(FL, ABF, color = Species)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ Species, scales = "free", ncol = 1)



###############
### Fit GAM ###
###############

dat2 <- dat |>
  nest(.by = Species)

# Fit models
dat3 <- dat2 |>
  mutate(fit = map(data, ~gam(ABF ~ s(FL, k = 5), data = .x, method = 'REML'))
         )

# Check goodness-of-fit
map(dat3$fit, gam.check)


# Calculate 1st order derivatives and find max (equivalent to root of 2nd derivative)
dat4 <- dat3 |>
  mutate(deriv = map(fit, ~derivatives(object = .x, term = "FL", order = 1,
                                       partial_match = TRUE))) |>
  mutate(inflect = map(deriv, ~{
    .x[which.max(.x$derivative), "data"]
    } |>
      as.numeric() |>
      round(2)
    ))
# Bull: none
# Blacktip: 103.5 cm
# Bonnethead: 79.3 cm



# Viz ABF vs FL w/ thresholds

inflect.df <- data.frame(Species = unique(dat$Species)[1:2] ,
                         FL = unlist(dat4$inflect)[1:2])

ggplot(data = dat, aes(FL, ABF, color = Species, fill = Species)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(formula = y ~ s(x, k = 5),
              method = "gam", linewidth = 0.7,
              method.args = list(family = gaussian(link = "identity")),
              alpha = 0.2
              ) +
  scale_color_met_d(palette_name = "Egypt", guide = "none") +
  scale_fill_met_d(palette_name = "Egypt", guide = "none") +
  geom_vline(data = inflect.df, aes(xintercept = FL), color = "black", linewidth = 1) +
  theme_bw(base_size = 16) +
  labs(x = "FL (cm)", y = "ABF (N)") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.margin=unit(rep(1, 4),'lines')) +
  facet_wrap(~ Species, scales = "free", ncol = 1)


# ggsave("Figures/Figure S2.tiff", width = 5, height = 7, units = "in", dpi = 400)
