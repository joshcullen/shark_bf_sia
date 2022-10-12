

# Function to calculate bilateral anterior bite force from digitized shark muscles
calc_bf <- function(file) {

  print(file)

# Load data
bite <- read.csv(file)
bite$Muscle[bite$Muscle == ""] <- NA

# Reformat data and create column for denoting muscle origin/insertion
bite <- bite %>%
  fill(Muscle) %>%
  mutate(loc = case_when(str_detect(string = Muscle, pattern = "-I$") ~ "Insertion",
                         str_detect(string = Muscle, pattern = "-O$") ~ "Origin",
                         TRUE ~ NA_character_),
         .after = Muscle) %>%
  mutate(Muscle = str_replace(string = Muscle, pattern = "-.", replacement = ""))


# Calculate mean position of each digitized location (of the 3 measurements per location)
coords.mean <- bite %>%
  group_by(Muscle, loc) %>%
  summarize(x = mean(X),
            y = mean(Y),
            z = mean(Z)) %>%
  ungroup()



#Muscle Distances (dm)
muscle.dist <- coords.mean %>%
  split(.$Muscle) %>%
  map(~euc_dist(start = .x[.x$loc == "Origin", c("x","y","z")],
                end = .x[.x$loc == "Insertion", c("x","y","z")])) %>%
  bind_rows(.id = "Muscle") %>%
  drop_na() %>%
  rename(dm = x)


#In-lever Distances (dil)
inlever.dist <- coords.mean %>%
  split(.$Muscle) %>%
  map(~euc_dist(start = coords.mean %>%
                  filter(Muscle == "JJ") %>%
                  dplyr::select(x, y, z),
                end = .x[.x$loc == "Insertion", c("x","y","z")])) %>%
  bind_rows(.id = "Muscle") %>%
  drop_na() %>%
  rename(dil = x)



#Muscle anatomical cross-sectional area
acsa <- bite %>%
  drop_na(A.CSA) %>%
  dplyr::select(Muscle, A.CSA) %>%
  mutate(fiber.type = case_when(Muscle == "QD4" ~ "Red",
                                TRUE ~ "White"))


#Individual Muscle Force Production (Po)--Unilateral
muscle.force <- acsa %>%
  mutate(force = case_when(fiber.type == "White" ~ A.CSA * 28.9,
                           fiber.type == "Red" ~ A.CSA * 14.2)) %>%
  arrange(Muscle)
#values for red and white muscle spec tension from Lou et al 2002

Ftot <- sum(muscle.force$force)  #total unilateral force



#Proportion of total tetanic force (Po) - used for in-lever calc
muscle.force <- muscle.force %>%
  mutate(force.prop = force / Ftot)


# Resolved in-lever (weighted) across all muscles
muscle.force <- muscle.force %>%
  left_join(., inlever.dist, by = "Muscle") %>%
  mutate(w.dil = force.prop * dil)


#Resolved jaw adductor in-lever
JAdil <- sum(muscle.force$w.dil)


#Anterior Bite Point Out-lever
ABPdol <- euc_dist(start = coords.mean %>%
                     filter(Muscle == "JJ") %>%
                     dplyr::select(x, y, z),
                   end = coords.mean %>%
                     filter(Muscle == "ABP") %>%
                     dplyr::select(x, y, z))

#Posterior Bite Point Out-lever
# PBPdol <- euc_dist(start = coords.mean %>%
#                      filter(Muscle == "JJ") %>%
#                      dplyr::select(x, y, z),
#                    end = coords.mean %>%
#                      filter(Muscle == "PBP") %>%
#                      dplyr::select(x, y, z))


#Anterior Mechanical Advantage (AMA)
AMA<- JAdil/ABPdol

#Posterior Mechanical Advantage (PMA)
# PMA<- JAdil/PBPdol


## Calculation of moments of individual muscles

#dummy axis of jaw joints
JJopp <- coords.mean[coords.mean$Muscle == "JJ", c('x','y','z')]
JJopp[2] <- JJopp[2] - 15  #right JJ y-coord (-15) chosen at random
rJJ <- (JJopp - coords.mean[coords.mean$Muscle == "JJ", c('x','y','z')]) / 15  #unit vector of JJ axis


# Muscle unit vectors (insertion coord-origin coord)/magnitude
muscle.uv <- coords.mean %>%
  drop_na(loc) %>%
  arrange(desc(loc)) %>%
  split(.$Muscle) %>%
  map(., ~{.x %>%
      dplyr::select(x, y, z) %>%
      apply(., 2, diff)}) %>%
  map2(.x = ., .y = muscle.dist$dm,
       ~{.x / .y}) %>%
  bind_rows(.id = "Muscle")

#Verify that unit vectors (root of sum of squared values equals 1)
apply(muscle.uv[,-1], 1, function(a) sqrt(sum(a^2)))


#creation of muscle force vectors (unit vector*force)
Fm <- muscle.uv %>%
  split(.$Muscle) %>%
  map2(.x = ., .y = muscle.force$force, ~{.x %>%
      mutate(across(x:z, function(a) a * .y))}) %>%
  bind_rows(.id = "Muscle")


# Calculating position vector of muscle insertion from JJ
rM <- coords.mean %>%
  filter(loc == "Insertion") %>%
  split(.$Muscle) %>%
  map(., ~{.x %>%
      dplyr::select(x, y, z)}) %>%
  map(., function(a) a - coords.mean[coords.mean$Muscle == "JJ", c("x","y","z")]) %>%
  bind_rows(.id = "Muscle")


# Create perpendicular vector to that of ABP/PBP-JJ-JJopp
rABP <- coords.mean %>%
  filter(Muscle %in% c('ABP','JJ')) %>%
  arrange(desc(Muscle)) %>%
  dplyr::select(x, y, z) %>%
  apply(., 2, diff) %>%
  as.matrix() %>%
  t()
# rPBP <- coords.mean %>%
#   filter(Muscle %in% c('PBP','JJ')) %>%
#   dplyr::select(x, y, z) %>%
#   apply(., 2, diff) %>%
#   as.matrix() %>%
#   t()




# Cross product of unit vectors results in perp vector to describe BF direction
rABF<- xprod(rABP, as.numeric(rJJ)) %>%
  t()  #transpose from col to row vector
uvABF <- rABF / vect_mag(rABF)

# rPBF<- xprod(rPBP, as.numeric(rJJ)) %>%
#   t()  #transpose from col to row vector
# uvPBF <- rPBF / vect_mag(rPBF)


#Proj(fv)=(fv dot uvABF)*uvABF for all muscle fvs

orthog.ABF <- apply(Fm[,-1], 1, function(a) a %*% t(uvABF))
biABF<- sum(orthog.ABF) * JAdil / ABPdol * 2 #bilateral ABF

# orthog.PBF <- apply(Fm[,-1], 1, function(a) a %*% t(uvPBF))
# biPBF<- sum(orthog.PBF) * JAdil / PBPdol * 2 #bilateral ABF




### Export results ###
id <- gsub(pattern = "Raw_data/Digitized Points/|_DigPts.csv", replacement = "", x = file)
species <- ifelse(grepl(pattern = "Cleu", x = id), "Bull",
                  ifelse(grepl(pattern = "BT", x = id), "Blacktip", "Bonnethead"))

out <- data.frame(SharkID = id,
                  Species = species,
                  ABF = as.numeric(biABF$x))

return(out)
}
