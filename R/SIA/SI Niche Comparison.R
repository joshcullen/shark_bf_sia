

#### Comparison of ecological niches based upon bulk carbon and nitrogen stable isotopes ####

library(tidyverse)
library(SIBER)
library(ggdist)
library(MetBrewer)


### Load data ###
sia <- read.csv("Raw_data/SIA Master.csv")

glimpse(sia)
summary(sia)
head(sia)
# sia$Sample.Yr<- factor(sia$Sample.Yr) #turn sampled year into a factor



### Clean and wrangle data ###

#Check for duplicates and remove record(s) for any duplicate IDs
ind <- which(duplicated(sia$SharkID) | is.na(sia$TL))
sia <- sia[-ind,]



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
siberDensityPlot(SEA.B, xticklabels = c('Cleu','Clim','Stib'),
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
  mutate(across(species, factor, level = c('Bull','Blacktip','Bonnethead')))


ggplot(SEA.B.df, aes(species, SEA)) +
  ggdist::stat_halfeye(aes(fill = species), adjust = 0.5, width = 0.6, .width = 0,
                       justification = -0.3, point_color = NA) +
  scale_fill_met_d(name = "Egypt") +
  geom_jitter(aes(color = species), width = .05, alpha = .05) +
  scale_color_met_d(name = "Egypt") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "transparent") +
  scale_y_continuous(breaks = seq(0, 14, by = 2)) +
  labs(x = "", y = expression("Standard Ellipse Area " ('\u2030' ^2) )) +
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




### Compare isotopic niche width ###

## Calculate probability that posterior niche width of one species is larger than another

# Bull vs blacktip
sum(SEA.B[,1] > SEA.B[,2]) / nrow(SEA.B)  #Bull shark has larger niche w/ 100% probability

# Bull vs bonnethead
sum(SEA.B[,1] > SEA.B[,3]) / nrow(SEA.B)  #Bull shark has larger niche w/ 100% probability

# Bonnethead vs blacktip
sum(SEA.B[,3] > SEA.B[,2]) / nrow(SEA.B)  #Bonnethead shark has larger niche w/ 100% probability





### Calculate Bayesian niche overlap ###

# Calc overlap of 95% ellipses
overlap.cleu.clim <- bayesianOverlap("1.Cleu", "1.BT", ellipses.posterior,
                                       draws = 10, p.interval = 0.95,
                                       n = 360)
overlap.cleu.clim <- overlap.cleu.clim %>%
  mutate(cleu_onto_clim = overlap / area2,
         clim_onto_cleu = overlap / area1)
summarize(overlap.cleu.clim,
          mean.cleu.clim = mean(cleu_onto_clim),
          sd.cleu.clim = sd(cleu_onto_clim),
          mean.clim.cleu = mean(clim_onto_cleu),
          sd.clim.cleu = sd(clim_onto_cleu))

overlap.cleu.stib <- bayesianOverlap("1.Cleu", "1.BH", ellipses.posterior,
                                     draws = 10, p.interval = 0.95,
                                     n = 360)
overlap.cleu.stib <- overlap.cleu.stib %>%
  mutate(cleu_onto_stib = overlap / area2,
         stib_onto_cleu = overlap / area1)
summarize(overlap.cleu.stib,
          mean.cleu.stib = mean(cleu_onto_stib),
          sd.cleu.stib = sd(cleu_onto_stib),
          mean.stib.cleu = mean(stib_onto_cleu),
          sd.stib.cleu = sd(stib_onto_cleu))

overlap.clim.stib <- bayesianOverlap("1.BT", "1.BH", ellipses.posterior,
                                     draws = 10, p.interval = 0.95,
                                     n = 360)
overlap.clim.stib <- overlap.clim.stib %>%
  mutate(clim_onto_stib = overlap / area2,
         stib_onto_clim = overlap / area1)
summarize(overlap.clim.stib,
          mean.clim.stib = mean(clim_onto_stib),
          sd.clim.stib = sd(clim_onto_stib),
          mean.stib.clim = mean(stib_onto_clim),
          sd.stib.clim = sd(stib_onto_clim))
