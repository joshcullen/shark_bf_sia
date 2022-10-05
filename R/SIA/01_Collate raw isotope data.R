

### Collate raw stable isotope data on Galveston Bay sharks ###

library(tidyverse)



t1 <- read.csv("Raw_data/Tray 1_7.5.16.csv")
t2 <- read.csv("Raw_data/Tray 2_7.18.16.csv")
t3 <- read.csv("Raw_data/Tray 1_2.20.18.csv")

SIA <- rbind(t1,t2,t3)
names(SIA)[1] <- "SharkID"

SIA$SharkID <- gsub("/", "-", SIA$SharkID)

#Check C:N
SIA$C.N <- SIA$wt.per.C/SIA$wt.per.N  #calc C:N for each sample

#do any samples have a C:N > 3.5?
SIA$C.N[SIA$C.N >= 3.5] #no, all samples < 3.5


# Create 'Species' column

SIA <- SIA %>%
  mutate(Species = case_when(str_detect(string = SharkID, pattern = "Cleu") ~ "Cleu",
                             str_detect(string = SharkID, pattern = "BT") ~ "BT",
                             str_detect(string = SharkID, pattern = "BH") ~ "BH"),
         .before = SharkID)





# Export collated master .csv file
write.csv(SIA, "Raw_data/SIA Master2.csv", row.names = F)
