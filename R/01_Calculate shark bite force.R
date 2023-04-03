
### Calculate theoretical bite force across shark species ###

library(tidyverse)
library(vroom)

source("R/helper_functions.R")
source("R/calc_bf.R")


# Load data files and metadata
files <- list.files("Raw_data/Digitized Points", full.names = TRUE)
metadata <- list.files("Raw_data", pattern = "metadata.csv$", full.names = TRUE) %>%
  vroom(., delim = ",")

# Calculate bite force
dat <- map(files, calc_bf) %>%
  bind_rows()

# Clean SharkID names
dat <- dat %>%
  mutate(SharkID = str_replace(SharkID, "_DigPts_redo.csv", "")) %>%
  mutate(SharkID = str_replace(SharkID, "_redo", "")) %>%
  mutate(SharkID = str_replace(SharkID, "_", " ")) %>%
  mutate(SharkID = str_replace_all(SharkID, "\\.", "-")) %>%
  mutate(SharkID = str_replace(SharkID, "Ara_15", "Ara15"))

metadata <- metadata %>%
  mutate(SharkID = str_replace(SharkID, "_", " "))



# Join metadata to bite force estimates
dat2 <- left_join(dat, metadata, by = "SharkID") %>%
  arrange(Species, TL)



# Viz BF over total length by species
ggplot(dat2, aes(TL, ABF, color = Species)) +
  geom_point() +
  geom_smooth(method = "gam") +
  theme_bw() +
  facet_wrap(~ Species, scales = "free", nrow = 3)



### Export results ###

write.csv(dat2, "Processed_data/Bite force results.csv", row.names = FALSE)
