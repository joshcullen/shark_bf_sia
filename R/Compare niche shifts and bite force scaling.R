
### Evaluate ontogenetic shifts in isotopic niche related to bite force scaling ###

library(tidyverse)
library(patchwork)
library(MetBrewer)
library(ggh4x)
library(mgcv)


#################
### Load data ###
#################

bf_dat <- read_csv("Processed_data/Bite force results.csv")
bf_dat <- bf_dat |>
  mutate(SharkID = case_when(str_detect(SharkID, "Cleu") ~ gsub(" ", "_", SharkID),
                             TRUE ~ SharkID))

sia_dat <- read_csv("Raw_data/SIA Master.csv")
sia_dat$Species <- factor(sia_dat$Species, level = c('Cleu','BT','BH'))
levels(sia_dat$Species) <- c('Bull','Blacktip','Bonnethead')


# Join datasets
dat <- full_join(bf_dat,
                 sia_dat[,c("SharkID", "Species", "FL", "d13C", "d15N")],
                 by = c("SharkID", "Species", "FL")) |>
  drop_na(FL)


# Define color palette for viz
pal <- met.brewer(palette_name = "Egypt", n = 3, return_hex = TRUE)



####################################################
### Assess ontogenetic niche shifts via isotopes ###
####################################################

# Fit models
dat_mod <- dat |>
  nest(.by = Species) |>
  mutate(d13C = map(data, ~gam(d13C ~ s(FL, k = 5), data = .x, method = 'REML')),
         d15N = map(data, ~gam(d15N ~ s(FL, k = 5), data = .x, method = 'REML')))

# Determine which relationships are significant
dat_mod2 <- dat_mod |>
  pivot_longer(cols = -c(Species, data), names_to = "mod", values_to = "fit") |>
  mutate(signif = map(fit, ~{summary(.x)$s.pv} |>
                        unlist() |>
                        round(4))) |>
  filter(signif < 0.05)

# Subset data for plotting GAMs
gam_smooth <- dat |>
  select(Species, FL, ABF, d13C, d15N) |>
  pivot_longer(cols = -c(Species, FL)) |>
  filter(Species == 'Blacktip' & name %in% c('d13C','d15N') |
           Species == 'Bonnethead' & name == 'd13C' |
           Species == 'Bull' & name == 'd13C') |>
  drop_na() |>
  mutate(Species = factor(Species, levels = c('Bull','Blacktip','Bonnethead')),
         name = case_when(name == 'd13C' ~ 'paste(delta^{13}, "C (\u2030)")',
                          name == 'd15N' ~ 'paste(delta^{15}, "N (\u2030)")',
                          name == 'ABF' ~ 'paste("ABF (N)")')) |>
  mutate(name = factor(name, levels = c('paste(delta^{13}, "C (\u2030)")',
                                        'paste(delta^{15}, "N (\u2030)")',
                                        'paste("ABF (N)")')))

levels(gam_smooth$Species) <- c("bold(Bull)", "bold(Blacktip)", "bold(Bonnethead)")




###############################
### Visualize relationships ###
###############################

# Reformat data
dat2 <- dat |>
  select(Species, FL, ABF, d13C, d15N) |>
  pivot_longer(cols = -c(Species, FL)) |>
  mutate(Species = factor(Species, levels = c('Bull','Blacktip','Bonnethead')),
         name = case_when(name == 'd13C' ~ 'paste(delta^{13}, "C (\u2030)")',
                          name == 'd15N' ~ 'paste(delta^{15}, "N (\u2030)")',
                          name == 'ABF' ~ 'paste("ABF (N)")')) |>
  mutate(name = factor(name, levels = c('paste(delta^{13}, "C (\u2030)")',
                                        'paste(delta^{15}, "N (\u2030)")',
                                        'paste("ABF (N)")')))

levels(dat2$Species) <- c("bold(Bull)", "bold(Blacktip)", "bold(Bonnethead)")




# Add color to species titles in plot
strip_cols <- strip_themed(
  text_x = elem_list_text(colour = pal),
  by_layer_x = FALSE
  )

ggplot() +
  geom_point(data = dat2, aes(FL, value, color = Species), size = 2, alpha = 0.7) +
  geom_point(data = dat2 |>
               filter(name == 'paste("ABF (N)")'),
             aes(FL, value), size = 2) +
  scale_color_met_d(palette_name = "Egypt", guide = "none") +
  geom_smooth(data = gam_smooth, mapping = aes(FL, value, color = Species, fill = Species),
              formula = y ~ s(x, k = 5),
              method = "gam", linewidth = 0.7,
              method.args = list(family = gaussian(link = "identity")),
              alpha = 0.2
  ) +
  scale_fill_met_d(palette_name = "Egypt", guide = "none") +
  labs(x = "FL (cm)", y = "") +
  theme_bw(base_size = 14) +
  facet_grid2(name ~ Species, scales = "free", independent = "y", switch = "y",
              labeller = label_parsed, strip = strip_cols) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14),
        strip.text.x = element_text(color = "red"),
        panel.grid = element_blank())

ggsave("Figures/Isotopes and BF over ontogeny.png", width = 10, height = 8, units = "in", dpi = 600)



