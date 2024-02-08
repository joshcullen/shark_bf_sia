

#### Comparison of ecological niches based upon bulk carbon and nitrogen stable isotopes ####

library(tidyverse)
library(SIBER)
library(ggdist)
library(MetBrewer)
library(patchwork)
library(tictoc)


### Load data ###
sia <- read.csv("Raw_data/SIA Master.csv")

glimpse(sia)
summary(sia)
head(sia)
sia$Species <- factor(sia$Species, level = c('Cleu','BT','BH'))
levels(sia$Species) <- c('Bull','Blacktip','Bonnethead')


### Clean and wrangle data ###

#Check for duplicates and remove record(s) for any duplicate IDs
ind <- which(duplicated(sia$SharkID) | is.na(sia$TL))
# sia <- sia[-ind,]



### Calculate summary niche metrics ###

# Mean and SD of d13C and d15N by spp
sia %>%
  group_by(Species) %>%
  summarize(mean.C = mean(d13C),
            sd.C = sd(d13C),
            mean.N = mean(d15N),
            sd.N = sd(d15N),
            n = n())




### Calculate isotopic niche width ###

# Prep object for SIBER
sia2 <- sia[,c('d13C','d15N','Species')] %>%
  rename(iso1 = d13C, iso2 = d15N, group = Species) %>%
  mutate(community = 1)

sia.siber <- createSiberObject(sia2)
sia.siber


plotSiberObject(sia.siber,
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'))



# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 3        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the
# means. Fitting is via the JAGS method.
set.seed(2022)
ellipses.posterior <- siberMVN(sia.siber, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group
SEA.B <- siberEllipses(ellipses.posterior)


# Viz comparison of isotopic niche areas
siberDensityPlot(SEA.B, xticklabels = c('Bull', 'Blacktip','Bonnethead'),
                 xlab = c("Species"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

SEA.B.df <- data.frame(SEA.B) %>%
  rename(Bull = X1, Blacktip = X2, Bonnethead = X3) %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(cols = -iter, names_to = "species", values_to = 'SEA') %>%
  mutate(across(species, \(x) factor(x, level = c('Bull','Blacktip','Bonnethead'))
                )
         )


p.nw <- ggplot(SEA.B.df, aes(species, SEA)) +
  ggdist::stat_halfeye(aes(fill = species), adjust = 0.5, width = 0.6, .width = 0,
                       justification = -0.3, point_color = NA) +
  scale_fill_met_d(palette_name = "Egypt") +
  geom_jitter(aes(color = species), width = .05, alpha = .05) +
  scale_color_met_d(palette_name = "Egypt") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "transparent") +
  scale_y_continuous(breaks = seq(0, 14, by = 2)) +
  labs(x = "", y = expression("Ellipse Area " ('\u2030' ^2) )) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 12))


#Calculate some credible intervals
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B),
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

SEA.B.credibles





### Plot subset Bayesian posterior ellipses ###

n.posts <- 50  #number of posterior draws
p.ell <- 0.95  #95% quantile ellipse


# a list to store the results
all_ellipses <- list()

# loop over groups
for (i in 1:length(ellipses.posterior)){
  print(i)

  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL

  for (j in 1:n.posts){

    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)

    # mean
    mu  <- ellipses.posterior[[i]][j,5:6]

    # ellipse points
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)


    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))

  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}

# Merge all species together
ellipse_df <- all_ellipses %>%
  bind_rows(.id = "species") %>%
  rename(d13C = x, d15N = y) %>%
  mutate(species = case_when(species == 1 ~ "Bull",
                             species == 2 ~ "Blacktip",
                             species == 3 ~ "Bonnethead")) %>%
  mutate(across(species, factor, level = c('Bull','Blacktip','Bonnethead')))



p.cn <- ggplot() +
  geom_point(data = sia %>%
               rename(species = Species), aes(d13C, d15N, color = species), size = 2, alpha = 0.7) +
  geom_path(data = ellipse_df %>%
              filter(species == "Bull"), aes(d13C, d15N, group = rep, color = species),
               alpha = 0.15) +
  geom_path(data = ellipse_df %>%
              filter(species == "Blacktip"), aes(d13C, d15N, group = rep, color = species),
            alpha = 0.15) +
  geom_path(data = ellipse_df %>%
              filter(species == "Bonnethead"), aes(d13C, d15N, group = rep, color = species),
            alpha = 0.15) +
  scale_color_met_d(palette_name = "Egypt") +
  labs(x = expression(paste(delta^{13}, "C (\u2030)")), y = expression(paste(delta^{15}, "N (\u2030)"))) +
  theme_bw() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid = element_blank()) +
  guides(color = "none")





## Create composite plot for SIA
p.cn / p.nw +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = 'top',
        # plot.tag.position = c(0.09, 1),
        plot.tag = element_text(size = 18, hjust = 0, vjust = 0, face = 'bold'))

# ggsave("Figures/SIA niche width.png", width = 6, height = 8, units = "in", dpi = 600)




### Compare isotopic niche width ###

## Calculate probability that posterior niche width of one species is larger than another

# Bull vs blacktip
sum(SEA.B[,3] > SEA.B[,1]) / nrow(SEA.B)  #Bull shark has larger niche w/ 100% probability

# Bull vs bonnethead
sum(SEA.B[,3] > SEA.B[,2]) / nrow(SEA.B)  #Bull shark has larger niche w/ 100% probability

# Bonnethead vs blacktip
sum(SEA.B[,2] > SEA.B[,1]) / nrow(SEA.B)  #Bonnethead shark has larger niche w/ 100% probability





### Calculate Bayesian niche overlap ###

# Calc overlap of 95% ellipses
tic()
overlap.cleu.clim <- bayesianOverlap("1.Bull", "1.Blacktip", ellipses.posterior,
                                       p.interval = 0.95, draws = 1000, n = 100)
toc()  #took 1 min for 1000 draws
BRRR::skrrrahh(sound = "khaled3")

overlap.cleu.clim <- overlap.cleu.clim %>%
  mutate(cleu_onto_clim = overlap / area2,
         clim_onto_cleu = overlap / area1)
summarize(overlap.cleu.clim,
          mean.cleu.clim = mean(cleu_onto_clim),
          sd.cleu.clim = sd(cleu_onto_clim),
          mean.clim.cleu = mean(clim_onto_cleu),
          sd.clim.cleu = sd(clim_onto_cleu))

tic()
overlap.cleu.stib <- bayesianOverlap("1.Bull", "1.Bonnethead", ellipses.posterior,
                                     p.interval = 0.95, draws = 1000, n = 100)
toc()  #took 1 min for 1000 draws
BRRR::skrrrahh(sound = "khaled3")

overlap.cleu.stib <- overlap.cleu.stib %>%
  mutate(cleu_onto_stib = overlap / area2,
         stib_onto_cleu = overlap / area1)
summarize(overlap.cleu.stib,
          mean.cleu.stib = mean(cleu_onto_stib),
          sd.cleu.stib = sd(cleu_onto_stib),
          mean.stib.cleu = mean(stib_onto_cleu),
          sd.stib.cleu = sd(stib_onto_cleu))

tic()
overlap.clim.stib <- bayesianOverlap("1.Blacktip", "1.Bonnethead", ellipses.posterior,
                                     p.interval = 0.95, draws = 1000, n = 100)
toc()  #took 1 min for 1000 draws
BRRR::skrrrahh(sound = "khaled3")

overlap.clim.stib <- overlap.clim.stib %>%
  mutate(clim_onto_stib = overlap / area2,
         stib_onto_clim = overlap / area1)
summarize(overlap.clim.stib,
          mean.clim.stib = mean(clim_onto_stib),
          sd.clim.stib = sd(clim_onto_stib),
          mean.stib.clim = mean(stib_onto_clim),
          sd.stib.clim = sd(stib_onto_clim))



# Create viz of overlap

cleu.over <- data.frame(From = "Bull",
                   To = rep(c("Blacktip", "Bonnethead"), each = nrow(overlap.cleu.clim)),
                   overlap = c(overlap.cleu.clim$cleu_onto_clim, overlap.cleu.stib$cleu_onto_stib))
clim.over <- data.frame(From = "Blacktip",
                        To = rep(c("Bull", "Bonnethead"), each = nrow(overlap.cleu.clim)),
                        overlap = c(overlap.cleu.clim$clim_onto_cleu, overlap.clim.stib$clim_onto_stib))
stib.over <- data.frame(From = "Bonnethead",
                        To = rep(c("Bull", "Blacktip"), each = nrow(overlap.cleu.clim)),
                        overlap = c(overlap.cleu.stib$stib_onto_cleu, overlap.clim.stib$stib_onto_clim))

over <- rbind(cleu.over, clim.over, stib.over) %>%
  mutate(across(From:To, factor, level = c('Bull','Blacktip','Bonnethead')))

med.over <- over %>%
  group_by(From, To) %>%
  summarize(value = median(overlap)) %>%
  ungroup()


ggplot(data = over, aes(overlap, fill = From)) +
  geom_density(aes(color = From), size = 0.75) +
  geom_vline(data = med.over, aes(xintercept = value), linewidth = 0.5) +
  scale_fill_met_d(palette_name = "Egypt") +
  scale_color_met_d(palette_name = "Egypt") +
  labs(x = "Overlap Probability", y = "Density") +
  scale_y_continuous(position = "right") +
  scale_x_continuous(limits = c(0, 1.04), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(fill = "none", color = "none") +
  facet_grid(From ~ To, scales = "free_y", switch = "y")

# ggsave("Figures/SIA niche overlap.png", width = 12, height = 8, units = "in", dpi = 600)
