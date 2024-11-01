

#### Comparison of ecological niches based upon bulk carbon and nitrogen stable isotopes ####

library(tidyverse)
library(nicheROVER)
library(ggdist)
library(MetBrewer)
library(patchwork)
library(tictoc)
library(ggh4x)


#################
### Load data ###
#################

sia <- read.csv("Raw_data/SIA Master.csv")

glimpse(sia)
summary(sia)
head(sia)

sia$Species <- factor(sia$Species, level = c('Cleu','BT','BH'))
levels(sia$Species) <- c('Bull','Blacktip','Bonnethead')


# Load metadata to determine seasons of samples
meta <- read_csv("Raw_data/Dock Sampling Inventory.csv")
meta2 <- meta %>%
  mutate(Date = gsub(pattern = "/", replacement = "-", x = Date),
         SharkID = case_when(str_detect(ID, "^Cleu") ~ ID,
                             TRUE ~ paste(ID, Date)),
         TL = as.numeric(TL)) %>%
  mutate(Date = as_date(Date, format = "%m-%d-%y"),
         month = month(Date),
         season = case_when(month %in% c(12, 1, 2) ~ 'Winter',
                            month %in% c(3:5) ~ 'Spring',
                            month %in% c(6:8) ~ 'Summer',
                            month %in% c(9:11) ~ 'Fall'),
         .after = Date)

# Add 'month' and 'season' to SIA data
sia2 <- left_join(sia,
                  meta2[,c("SharkID","TL","month","season")],
                  by = join_by("SharkID","TL"))



##############################
### Clean and wrangle data ###
##############################

#Check for duplicates and remove record(s) for any duplicate IDs
ind <- which(duplicated(sia2$SharkID) | is.na(sia2$TL))
sia2 <- sia2[-ind,]


#######################################
### Calculate summary niche metrics ###
#######################################

# Mean and SD of d13C and d15N by spp
sia2 %>%
  group_by(Species) %>%
  summarize(mean.C = mean(d13C),
            sd.C = sd(d13C),
            mean.N = mean(d15N),
            sd.N = sd(d15N),
            n = n())


# Mean and SD of d13C and d15N by spp and season
(sia_szn_n <- sia2 %>%
  group_by(Species, season) %>%
  summarize(mean.C = mean(d13C),
            sd.C = sd(d13C),
            mean.N = mean(d15N),
            sd.N = sd(d15N),
            n = n()))


#######################################################
### Explore size distribution by species and season ###
#######################################################

sia2 %>%
  mutate(season = factor(season, levels = c('Spring','Summer','Fall'))) %>%
  drop_na(season) %>%

ggplot(aes(TL)) +
  geom_density(aes(fill = Species)) +
  geom_rug() +
  scale_fill_met_d("Egypt", guide = "none") +
  labs(x = "TL (cm)", y = "Density") +
  theme_bw(base_size = 14) +
  theme(#axis.text = element_text(size=16),
        #axis.title = element_text(size=18),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold")) +
  facet_grid(season ~ Species, scales = "free_x")

ggsave("Figures/Figure S1.tiff", width = 8, height = 6, units = "in", dpi = 400)





######################################
### Calculate isotopic niche width ###
######################################

### For all samples per spp
# generate parameter draws from the "default" posteriors of each spp
sia3 <- sia2[,c('d13C','d15N','Species','season')]

set.seed(2024)
sia.par <- tapply(1:nrow(sia3), sia3$Species,
                  function(ii) niw.post(nsamples = 5e3, X = sia3[ii,1:2]))

niche.par.plot(sia.par, col = met.brewer('Egypt', n = 3), plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(sia.par), fill = met.brewer('Egypt', n = 3))



# Viz estimated niche regions
tmp <- tapply(1:nrow(sia3), sia3$Species,
                  function(ii) niw.post(nsamples = 50, X = sia3[ii,1:2]))
sia.data <- tapply(1:nrow(sia3), sia3$Species, function(ii) X = sia3[ii,1:2])

niche.plot(niche.par = tmp, niche.data = sia.data, pfrac = .05, alpha = 0.95,
           iso.names = expression(delta^{13}*C, delta^{15}*N),
           col = met.brewer('Egypt', n = 3), xlab = expression("Isotope Ratio (per mil)"))



### For each spp x season grouping
# generate parameter draws from the "default" posteriors of each spp x season
sia4 <- sia3 %>%
  drop_na(season)

set.seed(2024)
sia.par.szn <- tapply(1:nrow(sia4), paste(sia4$Species, sia4$season, sep = "_"),
                      function(ii) niw.post(nsamples = 5e3, X = sia4[ii,1:2]))



#####################################
### Create figures for manuscript ###
#####################################

### Plot subset of Bayesian posterior ellipses per species ###

n.post <- 50  #number of posterior draws

# define empty list to store resulting ellipses
all_ellipses <- vector("list", length = length(sia.par)) |>
  set_names(names(sia.par))

# generate ellipses from posterior
for (i in seq_along(sia.par)) {

  ell = NULL

  for (j in 1:n.post) {
    out <- nicheROVER::ellipse(mu = sia.par[[i]]$mu[j,], sia.par[[i]]$Sigma[,,j], alpha = 0.95, n = 100)
    out2 <- cbind(rep = j, out)

    ell <- rbind(ell, out2)
  }

  all_ellipses[[i]] <- data.frame(ell)
}


# Merge all species together
ellipse_df <- all_ellipses |>
  bind_rows(.id = "Species") |>
  rename(d13C = x, d15N = y) |>
  mutate(across(Species, \(x) factor(x, levels = c('Bull','Blacktip','Bonnethead'))))


# Create plot of ellipses
p.cn <- ggplot() +
  geom_point(data = sia, aes(d13C, d15N, color = Species), size = 2, alpha = 0.7) +
  geom_path(data = ellipse_df, aes(d13C, d15N, group = interaction(rep, Species), color = Species),
            alpha = 0.15) +
  scale_color_met_d("Egypt") +
  labs(x = expression(paste(delta^{13}, "C (\u2030)")), y = expression(paste(delta^{15}, "N (\u2030)"))) +
  theme_bw() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid = element_blank()) +
  guides(color = "none")



### Plot niche width ###

# posterior distribution of niche size by species
niche.width <- sapply(sia.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = 0.95)
})

# point estimate and standard error
rbind(est = colMeans(niche.width),
      se = apply(niche.width, 2, sd))


# convert niche width into long-format data.frame
niche.width.df <- data.frame(niche.width) |>
  mutate(tmp = 1:n()) |>
  pivot_longer(cols = -tmp, names_to = "Species", values_to = "niche_width") |>
  select(-tmp) |>
  mutate(Species = factor(Species, levels = c('Bull','Blacktip','Bonnethead')))


p.nw <- ggplot(niche.width.df, aes(Species, niche_width)) +
  ggdist::stat_halfeye(aes(fill = Species), adjust = 0.5, width = 0.6, .width = 0,
                       justification = -0.3, point_color = NA) +
  scale_fill_met_d("Egypt") +
  geom_jitter(aes(color = Species), width = .05, alpha = .05) +
  scale_color_met_d("Egypt") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "transparent") +
  # scale_y_continuous(breaks = seq(0, 14, by = 2)) +
  labs(x = "", y = expression("Ellipse Area " ('\u2030' ^2) )) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 12))





### Create composite plot ###

p.cn / p.nw +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = 'top',
        # plot.tag.position = c(0.09, 1),
        plot.tag = element_text(size = 18, hjust = 0, vjust = 0))

ggsave("Figures/Figure 4.tiff", width = 6, height = 8, units = "in", dpi = 400)







### Plot subset of Bayesian posterior ellipses per species x season ###

n.post <- 50  #number of posterior draws

# define empty list to store resulting ellipses
all_ellipses_szn <- vector("list", length = length(sia.par.szn)) |>
  set_names(names(sia.par.szn))

# generate ellipses from posterior
for (i in seq_along(sia.par.szn)) {

  ell = NULL

  for (j in 1:n.post) {
    out <- nicheROVER::ellipse(mu = sia.par.szn[[i]]$mu[j,],
                               sia.par.szn[[i]]$Sigma[,,j],
                               alpha = 0.95,
                               n = 100)
    out2 <- cbind(rep = j, out)

    ell <- rbind(ell, out2)
  }

  all_ellipses_szn[[i]] <- data.frame(ell)
}


# Merge all species together
ellipse_df_szn <- all_ellipses_szn |>
  bind_rows(.id = "group") |>
  separate(col = "group", into = c("Species",'season')) |>
  rename(d13C = x, d15N = y) |>
  mutate(Species = factor(Species, levels = c('Bull','Blacktip','Bonnethead')),
         season = factor(season, levels = c('Spring','Summer','Fall')))


# Create plot of ellipses
ggplot() +
  geom_point(data = sia4, aes(d13C, d15N, color = Species), size = 2, alpha = 0.7) +
  geom_path(data = ellipse_df_szn, aes(d13C, d15N, group = interaction(rep, Species),
                                       color = Species), alpha = 0.15) +
  scale_color_met_d("Egypt") +
  labs(x = expression(paste(delta^{13}, "C (\u2030)")), y = expression(paste(delta^{15}, "N (\u2030)"))) +
  coord_equal() +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = c(0.9,0.1),
        legend.justification = c(1,0)) +
  facet_wrap(~ season, nrow = 2)

ggsave("Figures/Figure S3.tiff", width = 6, height = 5, units = "in", dpi = 400)



### Plot niche width ###

# posterior distribution of niche size by species
niche.width.szn <- sapply(sia.par.szn, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = 0.95)
})

# point estimate and standard error
rbind(est = colMeans(niche.width.szn),
      se = apply(niche.width.szn, 2, sd))


# convert niche width into long-format data.frame
niche.width.szn.df <- data.frame(niche.width.szn) |>
  mutate(tmp = 1:n()) |>
  pivot_longer(cols = -tmp, names_to = "group", values_to = "niche_width") |>
  separate(col = "group", into = c("Species",'season')) |>
  select(-tmp) |>
  mutate(Species = factor(Species, levels = c('Bull','Blacktip','Bonnethead')),
         season = factor(season, levels = c('Spring','Summer','Fall')))


# create labels for N by spp x season
sia_szn_labs <- sia_szn_n |>
  drop_na(season) |>
  select(Species, season, n) |>
  mutate(y = -2)


ggplot(niche.width.szn.df, aes(season, niche_width)) +
  ggdist::stat_halfeye(aes(fill = Species), adjust = 0.5, width = 0.6, .width = 0,
                       justification = -0.3, point_color = NA) +
  scale_fill_met_d("Egypt", guide = "none") +
  geom_jitter(aes(color = Species), width = .05, alpha = .05) +
  scale_color_met_d("Egypt", guide = "none") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "transparent") +
  geom_text(data = sia_szn_labs,
            aes(season, y, label = paste0("n=",n)), size = 5) +
  labs(x = "", y = expression("Ellipse Area " ('\u2030' ^2) )) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold")) +
  facet_wrap(~ Species, nrow = 3, scales = "free_y")

ggsave("Figures/Figure S4.tiff", width = 6, height = 8, units = "in", dpi = 400)








####################################
### Compare isotopic niche width ###
####################################

## Calculate probability that posterior niche width of one species is larger than another

# Bull vs blacktip
apply(niche.width[,1:2], 1, which.max) |> table()  #Bull shark has larger niche w/ 100% probability

# Bull vs bonnethead
apply(niche.width[,c(1,3)], 1, which.max) |> table()  #Bull shark has larger niche w/ 100% probability

# Bonnethead vs blacktip
apply(niche.width[,2:3], 1, which.max) |> table()  #Bonnethead shark has larger niche w/ 100% probability




########################################
### Calculate Bayesian niche overlap ###
########################################

# Calc overlap of 95% ellipses
over_stat <- overlap(sia.par, nreps = 5000, nprob = 1000, alpha = 0.95)

over_stat_df <- over_stat |>
  as_tibble(rownames = "species_a") |>
  mutate(id = 1:n(),
         species_a = factor(species_a,
                            levels = c('Bull','Blacktip','Bonnethead'))) |>
  pivot_longer(cols = -c(id, species_a),
               names_to = "species_b",
               values_to = "over")  |>
  separate(species_b, into = c("species_c", "sample_number"),
           sep = "\\.") |>
  select(-id) |>
  rename(species_b = species_c) |>
  mutate(species_b =  factor(species_b,
                             level = c('Bull','Blacktip','Bonnethead')))


over_summ <- over_stat_df |>
  group_by(species_a, species_b) |>
  summarize(mean_over = round(mean(over), digits = 2),
            quant_2.5 = round(quantile(over, probs = 0.025, na.rm = TRUE), digits = 2),
            quant_97.5 = round(quantile(over, probs = 0.975, na.rm = TRUE), digits = 2)) |>
  ungroup() |>
  pivot_longer(cols = -c(species_a, species_b, mean_over),
               names_to = "percentage",
               values_to = "quant_val") |>
  mutate(
    percentage = as.numeric(str_remove(percentage, "quant_"))
  )


# Create viz of overlap
ggplot(data = over_stat_df, aes(over, fill = species_b)) +
  geom_density(aes(color = species_b), linewidth = 0.75) +
  geom_vline(data = over_summ, aes(xintercept = mean_over), linewidth = 0.75) +
  geom_vline(data = over_summ, aes(xintercept = quant_val),
             colour = "black", linewidth = 0.5, linetype = 6) +
  scale_fill_met_d(name = "Egypt") +
  scale_color_met_d(name = "Egypt") +
  labs(x = "Overlap Probability", y = "Density") +
  scale_x_continuous(limits = c(0, 1.04), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(fill = "none", color = "none") +
  ggh4x::facet_grid2(species_b ~ species_a,
                     independent = "y",
                     scales = "free_y")

ggsave("Figures/Figure 5.tiff", width = 12, height = 8, units = "in", dpi = 300)
