
### Evaluate ontogenetic shifts in isotopic niche related to bite force scaling ###

library(tidyverse)
library(patchwork)
library(MetBrewer)
library(ggh4x)
library(mgcv)
library(slider)



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
                 by = c("SharkID","Species","FL")) |>
  drop_na(FL)


# Define color palette for viz
pal <- met.brewer("Egypt", n = 3)



####################################################
### Assess ontogenetic niche shifts via isotopes ###
####################################################

# Fit models
dat_mod <- dat |>
  nest(.by = Species) |>
  mutate(d13C = map(data, ~gam(d13C ~ s(FL, k = 5), data = .x, method = 'REML')),
         d15N = map(data, ~gam(d15N ~ s(FL, k = 5), data = .x, method = 'REML')))

# Check goodness-of-fit
map(dat_mod$d13C, gam.check)
map(dat_mod$d15N, gam.check)

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
  # filter(Species == 'Blacktip' & name %in% c('d13C','d15N') |
  #          Species == 'Bonnethead' & name == 'd13C' |
  #          Species == 'Bull' & name == 'd13C') |>
  inner_join(dat_mod2[,c('Species','mod')],
             by = join_by("Species", "name" == "mod")) |>
  drop_na() |>
  mutate(Species = factor(Species, levels = c('Bull','Blacktip','Bonnethead')),
         name = case_when(name == 'd13C' ~ 'paste(delta^{13}, "C (\u2030)")',
                          name == 'd15N' ~ 'paste(delta^{15}, "N (\u2030)")',
                          name == 'ABF' ~ 'paste("ABF (N)")')) |>
  mutate(name = factor(name, levels = c('paste(delta^{13}, "C (\u2030)")',
                                        'paste(delta^{15}, "N (\u2030)")',
                                        'paste("ABF (N)")')))

levels(gam_smooth$Species) <- c("bold(Bull)", "bold(Blacktip)", "bold(Bonnethead)")




### Visualize relationships ###

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


# Store R^2 results for showing w/ signif model fits
r2_labs <- dat_mod2 |>
  # Create new columns for species-specific metrics
  mutate(R2 = map(fit, ~summary(.x)$r.sq),
         FL = map(data, ~{max(.x$FL) * 0.9}),
         value = case_when(mod == 'd13C' ~ map(data, ~{min(.x$d13C, na.rm = TRUE)}),
                           mod == 'd15N' ~ map(data, ~{min(.x$d15N, na.rm = TRUE)})),
         nudge = case_when(mod == 'd13C' ~ map(data, ~{abs(diff(range(.x$d13C, na.rm = TRUE))) / 24}),
                           mod == 'd15N' ~ map(data, ~{abs(diff(range(.x$d15N, na.rm = TRUE))) / 24}))) |>
  # Modify values as vector instead of list and add styling
  mutate(R2 = round(unlist(R2), 2),
         FL = unlist(FL),
         value = unlist(value),
         nudge = unlist(nudge),
         name = case_when(mod == 'd13C' ~ 'paste(delta^{13}, "C (\u2030)")',
                          mod == 'd15N' ~ 'paste(delta^{15}, "N (\u2030)")'),
         Species = paste0("bold(", Species, ")")) |>
  # Convert 'Species' and 'name' to factors
  mutate(Species = factor(Species, levels = levels(dat2$Species)),
         name = factor(name, levels = levels(dat2$name)),
         value = value + nudge) |>
  # Select only relevant columns
  select(Species, mod, FL, value, R2, name)



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
  scale_color_met_d("Egypt") +
  geom_smooth(data = gam_smooth, mapping = aes(FL, value, color = Species, fill = Species),
              formula = y ~ s(x, k = 5),
              method = "gam", linewidth = 0.7,
              method.args = list(family = gaussian(link = "identity")),
              alpha = 0.2
  ) +
  scale_fill_met_d("Egypt") +
  geom_text(data = r2_labs, aes(x = FL, y = value,
                                label = paste("italic(R)^2==", round(R2, 2))),
            parse = TRUE) +
  labs(x = "FL (cm)", y = "") +
  theme_bw(base_size = 14) +
  facet_grid2(name ~ Species, scales = "free", independent = "y", switch = "y",
              labeller = label_parsed, strip = strip_cols) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14),
        strip.text.x = element_text(color = "red"),
        panel.grid = element_blank(),
        legend.position = "none")

ggsave("Figures/Figure 3.tiff", width = 10, height = 8, units = "in", dpi = 300)






#######################################################
### Directly compare bite force and stable isotopes ###
#######################################################

### Calculate correlation using sliding window ###

dat4 <- dat |>
  group_by(Species) |>
  arrange(Species, FL) |>
  mutate(corr_d15N = slide_index2_dbl(.x = ABF,
                                 .y = d15N,
                                 .i = FL,
                                 .f = ~cor(.x, .y, use = "na.or.complete"),
                                 .before = 5,  #up to 5 cm FL less
                                 .after = 5),  #up to 5 cm FL more
         corr_d13C = slide_index2_dbl(.x = ABF,
                                .y = d13C,
                                .i = FL,
                                .f = ~cor(.x, .y, use = "na.or.complete"),
                                .before = 5,  #up to 5 cm FL less
                                .after = 5)  #up to 5 cm FL more
         ) |>
  ungroup()



# Viz relationships
dat4 |>
  select(Species, FL, corr_d15N, corr_d13C) |>
  pivot_longer(cols = c(corr_d13C, corr_d15N)) |>
  mutate(Species = factor(Species, levels = c("Bull", "Blacktip", "Bonnethead"))) |>

ggplot() +
  geom_line(aes(FL, value, color = Species), linewidth = 1.25) +
  geom_hline(yintercept = 0) +
  scale_color_met_d("Egypt", guide = "none") +
  labs(x = "FL (cm)", y = "Correlation", title = "Correlation between bite force and stable isotopes") +
  theme_bw(base_size = 14) +
  facet_grid(name ~ Species, scales = "free_x")
