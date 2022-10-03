
### Analyze bulk carbon and nitrogen stable isotopes for 3 shark spp ###


library(tidyverse)
# library(cowplot)

SIA <- read.csv("Raw_data/SIA Master.csv")

glimpse(SIA)
summary(SIA)
# SIA$Sample.Yr<- factor(SIA$Sample.Yr) #turn sampled year into a factor

#Check for duplicates and remove record(s) for any duplicate IDs
ind <- which(duplicated(SIA$SharkID) | is.na(SIA$TL))
SIA <- SIA[-ind,]


#To determine general metrics of sample sizes (for ABF and/or SIA); have a DF ready for ABF of all species together (ex. bfdata)


setwd("~/Documents/Shark Project/BF Data")

cleu.bfdata<- read.csv("Cleu BF Values.csv", header=TRUE, sep=",")[,1:23]
cleu.bfdata$Species<- rep("Cleu", nrow(cleu.bfdata))
clim.bfdata<- read.csv("BT BF Values.csv", header=TRUE, sep=",")[,1:23]
clim.bfdata$Species<- rep("Clim", nrow(clim.bfdata))
stib.bfdata<- read.csv("BH BF Values.csv", header=TRUE, sep=",")[,1:23]
stib.bfdata$Species<- rep("Stib", nrow(stib.bfdata))

bfdata<- rbind(cleu.bfdata, clim.bfdata, stib.bfdata)

BF_SIA.masterDF<- full_join(bfdata, SIA, by = "SharkID") #combine all individuals with some form of data being used

#Make full column of Species list
for(i in 1:nrow(BF_SIA.masterDF)) {
  if (is.na(BF_SIA.masterDF$Species.x[i])) {
    BF_SIA.masterDF$Species.x[i]<- paste0(BF_SIA.masterDF$Species.y[i])
  }
  else {
    BF_SIA.masterDF$Species.x[i]<- paste0(BF_SIA.masterDF$Species.x[i])
  }
}


#Make full column of Sex list
for(i in 1:nrow(BF_SIA.masterDF)) {
  if (is.na(BF_SIA.masterDF$Sex.x[i])) {
    BF_SIA.masterDF$Sex.x[i]<- paste0(BF_SIA.masterDF$Sex.y[i])
  }
  else {
    BF_SIA.masterDF$Sex.x[i]<- paste0(BF_SIA.masterDF$Sex.x[i])
  }
}


#For summarizing TL, FL, and PCL of all individuals

for(i in 1:nrow(BF_SIA.masterDF)) {
  if (is.na(BF_SIA.masterDF$TL.x[i])) {
    BF_SIA.masterDF$TL.x[i]<- paste0(BF_SIA.masterDF$TL.y[i])
  }
  else {
    BF_SIA.masterDF$TL.x[i]<- paste0(BF_SIA.masterDF$TL.x[i])
  }
}


for(i in 1:nrow(BF_SIA.masterDF)) {
  if (is.na(BF_SIA.masterDF$FL.x[i])) {
    BF_SIA.masterDF$FL.x[i]<- paste0(BF_SIA.masterDF$FL.y[i])
  }
  else {
    BF_SIA.masterDF$FL.x[i]<- paste0(BF_SIA.masterDF$FL.x[i])
  }
}


for(i in 1:nrow(BF_SIA.masterDF)) {
  if (is.na(BF_SIA.masterDF$PCL.x[i])) {
    BF_SIA.masterDF$PCL.x[i]<- paste0(BF_SIA.masterDF$PCL.y[i])
  }
  else {
    BF_SIA.masterDF$PCL.x[i]<- paste0(BF_SIA.masterDF$PCL.x[i])
  }
}

BF_SIA.masterDF$TL.x<- as.numeric(BF_SIA.masterDF$TL.x) #make col numeric
BF_SIA.masterDF$FL.x<- as.numeric(BF_SIA.masterDF$FL.x) #make col numeric
BF_SIA.masterDF$PCL.x<- as.numeric(BF_SIA.masterDF$PCL.x) #make col numeric


#Replace "BT" and "BH" with "Clim" and "Stib" for consistency of factors
BF_SIA.masterDF$Species.x<- gsub("BT", "Clim", BF_SIA.masterDF$Species.x)
BF_SIA.masterDF$Species.x<- gsub("BH", "Stib", BF_SIA.masterDF$Species.x)




#Summary stats

BF_SIA.masterDF %>% group_by(Species.x) %>% summarize(n = n()) # total N
BF_SIA.masterDF %>% filter(!is.na(ABF)) %>% group_by(Species.x) %>% summarize(n = n()) #BF n
BF_SIA.masterDF %>% filter(!is.na(d13C)) %>% group_by(Species.x) %>% summarize(n = n()) #SIA n

BF_SIA.masterDF %>% group_by(Species.x) %>% summarise(mean = mean(TL.x),
                                                      min = min(TL.x), max = max(TL.x)) #TL
BF_SIA.masterDF %>% group_by(Species.x) %>% summarise(mean = mean(FL.x),
                                                      min = min(FL.x), max = max(FL.x)) #FL

summary(BF_SIA.masterDF[BF_SIA.masterDF$Species.x == "Cleu",3]) #Cleu sex ratio
summary(BF_SIA.masterDF[BF_SIA.masterDF$Species.x == "Clim",3]) #Clim sex ratio
summary(BF_SIA.masterDF[BF_SIA.masterDF$Species.x == "Stib",3]) #Stib sex ratio



#Test for differences in d13C/d15N by year w/in each species
library(nlme) #for gls anova

ggplot(SIA, aes(x = Sample.Yr, y = d13C)) + geom_boxplot() + facet_wrap(~Species) #viz
ggplot(SIA, aes(x = Sample.Yr, y = d15N)) + geom_boxplot() + facet_wrap(~Species)

#Cleu
summary(SIA[SIA$Species == "Cleu",]$Sample.Yr) #use GLS Anova due to unequal n
shapiro.test(SIA[SIA$Species == "Cleu",]$d13C) #normal
shapiro.test(SIA[SIA$Species == "Cleu",]$d15N) #normal

Cleu.d13C.yr<- gls(d13C ~ Sample.Yr, data = SIA[SIA$Species == "Cleu",], weights = varIdent(form = ~ 1 | Sample.Yr), na.action = "na.omit")
anova(Cleu.d13C.yr, type = "marginal") #NS
summary(Cleu.d13C.yr)
plot(Cleu.d13C.yr) #check for HOV
qqnorm(Cleu.d13C.yr, abline = c(0,1)) #check normality


Cleu.d15N.yr<- gls(d15N ~ Sample.Yr, data = SIA[SIA$Species == "Cleu",], weights = varIdent(form = ~ 1 | Sample.Yr), na.action = "na.omit")
anova(Cleu.d15N.yr, type = "marginal") #NS
summary(Cleu.d15N.yr)
plot(Cleu.d15N.yr) #check for HOV
qqnorm(Cleu.d15N.yr, abline = c(0,1)) #check normality




#BT
summary(SIA[SIA$Species == "BT",]$Sample.Yr) #use GLS Anova due to unequal n
shapiro.test(SIA[SIA$Species == "BT",]$d13C) #non-normal
shapiro.test(SIA[SIA$Species == "BT",]$d15N) #normal

BT.d13C.yr<- gls(d13C ~ Sample.Yr, data = SIA[SIA$Species == "BT",], weights = varIdent(form = ~ 1 | Sample.Yr), na.action = "na.omit")
anova(BT.d13C.yr, type = "marginal") #NS
summary(BT.d13C.yr)
plot(BT.d13C.yr) #check for HOV
qqnorm(BT.d13C.yr, abline = c(0,1)) #check normality; close enough


BT.d15N.yr<- gls(d15N ~ Sample.Yr, data = SIA[SIA$Species == "BT",], weights = varIdent(form = ~ 1 | Sample.Yr), na.action = "na.omit")
anova(BT.d15N.yr, type = "marginal") #NS
summary(BT.d15N.yr)
plot(BT.d15N.yr) #check for HOV
qqnorm(BT.d15N.yr, abline = c(0,1)) #check normality






#BH
summary(SIA[SIA$Species == "BH",]$Sample.Yr) #use GLS Anova due to unequal n
shapiro.test(SIA[SIA$Species == "BH",]$d13C) #non-normal
shapiro.test(SIA[SIA$Species == "BH",]$d15N) #normal

BH.d13C.yr<- gls(d13C ~ Sample.Yr, data = SIA[SIA$Species == "BH",], weights = varIdent(form = ~ 1 | Sample.Yr), na.action = "na.omit")
anova(BH.d13C.yr, type = "marginal") #NS
summary(BH.d13C.yr)
plot(BH.d13C.yr) #check for HOV
qqnorm(BH.d13C.yr, abline = c(0,1)) #check normality; close enough


BH.d15N.yr<- gls(d15N ~ Sample.Yr, data = SIA[SIA$Species == "BH",], weights = varIdent(form = ~ 1 | Sample.Yr), na.action = "na.omit")
anova(BH.d15N.yr, type = "marginal") #NS
summary(BH.d15N.yr)
plot(BH.d15N.yr) #check for HOV
qqnorm(BH.d15N.yr, abline = c(0,1)) #check normality

######### No differences among years; pool all samples ##############





#Make biplots and look at trends over FL

#BH
ggplot(SIA[SIA$Species=="BH",], aes(d13C,d15N,color=Age,shape=Sex)) + geom_point(size=3)
ggplot(SIA[SIA$Species=="BH",], aes(FL,d15N,color=Age,shape=Sex)) + geom_point(size=3)
ggplot(SIA[SIA$Species=="BH",], aes(FL,d13C,color=Age,shape=Sex)) + geom_point(size=3)

#BT
ggplot(SIA[SIA$Species=="BT",], aes(d13C,d15N,color=Age,shape=Sex)) + geom_point(size=3)
ggplot(SIA[SIA$Species=="BT",], aes(FL,d15N,color=Age,shape=Sex)) + geom_point(size=3)
ggplot(SIA[SIA$Species=="BT",], aes(FL,d13C,color=Age,shape=Sex)) + geom_point(size=3)

#Cleu; only 2 females, so won't look at Sex diffs
ggplot(SIA[SIA$Species=="Cleu",], aes(d13C,d15N,color=Age)) + geom_point(size=3)
ggplot(SIA[SIA$Species=="Cleu",], aes(FL,d15N,color=Age)) + geom_point(size=3)
ggplot(SIA[SIA$Species=="Cleu",], aes(FL,d13C,color=Age)) + geom_point(size=3)

#Look at all 3 species together
ggplot(SIA, aes(d13C,d15N,color=Species)) + geom_point(size=3) + theme_bw()
ggplot(SIA, aes(FL,d15N,color=Species)) + geom_point(size=3) + geom_smooth(method = "lm", se=F) + theme_bw()
ggplot(SIA, aes(FL,d13C,color=Species)) + geom_point(size=3) + geom_smooth(method = "lm", se=F) + theme_bw()
ggplot(SIA, aes(ABF,d15N,color=Species)) + geom_point(size=3)



##Plotting with ABF and d15N dual axes

#Import BF dataset
setwd("~/Documents/Shark Project/BF Data")

cleu.bfdata<- read.csv("Cleu BF Values.csv", header=TRUE, sep=",")
clim.bfdata<- read.csv("BT BF Values.csv", header=TRUE, sep=",")
stib.bfdata<- read.csv("BH BF Values.csv", header=TRUE, sep=",")


#BH
par(mar=c(5,5,2,5))
with(stib.bfdata, plot(FL,ABF, pch=16, cex=1.5, col='goldenrod', xlab = "FL (cm)", ylab = "ABF (N)"))
par(new=T)
with(SIA[SIA$Species=="BH",], plot(FL,d15N, pch=16, cex=1, col='black', axes=F, xlab = NA, ylab=NA))
axis(side=4)
mtext(side=4, line=3, expression(δ^{15}*"N (‰)"))
legend("topleft", legend = c("ABF",expression(δ^{15}*"N")), pch = c(16,16), col = c("goldenrod","black"), bty = 'n')


#BT
par(mar=c(5,5,2,5))
with(clim.bfdata, plot(FL,ABF, pch=16, cex=1.5, col='indianred4', main = "BT Bite Force and d15N"))
par(new=T)
with(SIA[SIA$Species=="BT",], plot(FL,d15N, pch=16, cex=1, col='blue', axes=F, xlab = NA, ylab=NA))
axis(side=4)
mtext(side=4, line=3, "d15N")
legend("topleft", legend = c("ABF","d15N"), pch = c(16,16), col = c("indianred4","blue"), bty = 'n')


#Cleu
par(mar=c(5,5,2,5))
with(cleu.bfdata, plot(FL,ABF, pch=16, cex=1.5, col='darkblue', main = "Cleu Bite Force and d15N"))
par(new=T)
with(SIA[SIA$Species=="Cleu",], plot(FL,d15N, pch=16, cex=1, col='red', axes=F, xlab = NA, ylab=NA))
axis(side=4)
mtext(side=4, line=3, "d15N")
legend("topleft", legend = c("ABF","d15N"), pch = c(16,16), col = c("darkblue","red"), bty = 'n')










##Plotting with ABF and d13C dual axes


#BH
par(mar=c(5,5,2,5))
with(stib.bfdata, plot(FL,ABF, pch=16, cex=1.5, col='goldenrod', main = "BH Bite Force and d13C"))
par(new=T)
with(SIA[SIA$Species=="BH",], plot(FL,d13C, pch=16, cex=1, col='black', axes=F, xlab = NA, ylab=NA))
axis(side=4)
mtext(side=4, line=3, "d13C")
legend("topleft", legend = c("ABF","d13C"), pch = c(16,16), col = c("goldenrod","black"), bty = 'n')


#BT
par(mar=c(5,5,2,5))
with(clim.bfdata, plot(FL,ABF, pch=16, cex=1.5, col='indianred4', main = "BT Bite Force and d13C"))
par(new=T)
with(SIA[SIA$Species=="BT",], plot(FL,d13C, pch=16, cex=1, col='blue', axes=F, xlab = NA, ylab=NA))
axis(side=4)
mtext(side=4, line=3, "d13C")
legend("topleft", legend = c("ABF","d13C"), pch = c(16,16), col = c("indianred4","blue"), bty = 'n')


#Cleu
par(mar=c(5,5,2,5))
with(cleu.bfdata, plot(FL,ABF, pch=16, cex=1.5, col='darkblue', main = "Cleu Bite Force and d13C"))
par(new=T)
with(SIA[SIA$Species=="Cleu",], plot(FL,d13C, pch=16, cex=1, col='red', axes=F, xlab = NA, ylab=NA))
axis(side=4)
mtext(side=4, line=3, "d13C")
legend("bottomright", legend = c("ABF","d13C"), pch = c(16,16), col = c("darkblue","red"), bty = 'n')






### For plot of everything ###

#Re-order species for plot
SIA$Species<- factor(SIA$Species, levels = c("Cleu", "BT", "BH"))

ggplot(SIA, aes(x = d13C, y = d15N, color=Species)) + geom_point(alpha = 0.65, size=6) + theme_bw() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), legend.position = c(0.84,0.8), legend.background = element_blank(), legend.title = element_blank(), legend.key.size = unit(2.5, 'lines')) + labs(x = expression(δ^{13}*"C (‰)"), y = expression(δ^{15}*"N (‰)")) + scale_color_manual(values = c('cornflowerblue','indianred','darkseagreen'), labels = c("C. leucas", "C. limbatus", "S. tiburo")) + guides(color = guide_legend(override.aes = list(size=4), label.theme = element_text(face = "italic", size = 14)))

ggsave("SIA plot_test.png", width = 9, height = 9, units = "in", dpi = 600)









#### Single-tissue analysis of variation in d13C/d15N over FL by analyzing resids ####

##Cleu; only 2 females, don't need to test for sex differences

SIA %>% filter(Species == "Cleu") %>% select(d13C) %>% as.matrix() %>% shapiro.test() #normal
SIA %>% filter(Species == "Cleu") %>% select(d15N) %>% as.matrix() %>% shapiro.test() #normal


#d13C; linear regression fits better, but no trend variation over ontogeny
ggplot(SIA[SIA$Species=="Cleu",], aes(x=FL, y=d13C)) + geom_point(size = 3, color = "cornflowerblue") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

Cleu.d13C.FL.lm<- lm(d13C ~ FL, data = SIA[SIA$Species=="Cleu",])
Cleu.d13C.FL.poly<- lm(d13C ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="Cleu",])
summary(Cleu.d13C.FL.lm) # F=0.1275, R2=0.01, NS
summary(Cleu.d13C.FL.poly) # F=1.679, R2=0.13, NS
anova(Cleu.d13C.FL.lm, Cleu.d13C.FL.poly) #stick with poly


Cleu.d13C.FL.polyResids<- abs(residuals(Cleu.d13C.FL.poly)) #want to analyze magnitude (abs val) of resids

ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=FL, y=Cleu.d13C.FL.polyResids)) + geom_point(size = 3, color = "cornflowerblue") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

Cleu.d13Cresids.FL.lm<- lm(Cleu.d13C.FL.polyResids ~ FL, data = SIA[SIA$Species=="Cleu",])
Cleu.d13Cresids.FL.poly<- lm(Cleu.d13C.FL.polyResids ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="Cleu",])
summary(Cleu.d13Cresids.FL.lm) #NS; F=0.932, R2=0.04
summary(Cleu.d13Cresids.FL.poly) #NS; F=0.565, R2=0.05
anova(Cleu.d13Cresids.FL.lm, Cleu.d13Cresids.FL.poly) #lm is better


#d15N; quadratic regression fits better, which shows variance greatest at juvenile size class
ggplot(SIA[SIA$Species=="Cleu",], aes(x=FL, y=d15N)) + geom_point(size = 3, color = "cornflowerblue") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

Cleu.d15N.FL.lm<- lm(d15N ~ FL, data = SIA[SIA$Species=="Cleu",])
Cleu.d15N.FL.poly<- lm(d15N ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="Cleu",])
summary(Cleu.d15N.FL.lm) # F=1.911, R2=0.08, NS
summary(Cleu.d15N.FL.poly) # F=3.33, R2=0.23, NS
anova(Cleu.d15N.FL.lm, Cleu.d15N.FL.poly) #poly is better


Cleu.d15N.FL.polyResids<- abs(residuals(Cleu.d15N.FL.poly)) #want to analyze magnitude (abs val) of resids

ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=FL, y=Cleu.d15N.FL.polyResids)) + geom_point(size = 3, color = "cornflowerblue") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

Cleu.d15Nresids.FL.lm<- lm(Cleu.d15N.FL.polyResids ~ FL, data = SIA[SIA$Species=="Cleu",])
Cleu.d15Nresids.FL.poly<- lm(Cleu.d15N.FL.polyResids ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="Cleu",])
summary(Cleu.d15Nresids.FL.lm) #NS; F=0.0225, R2=0.001
summary(Cleu.d15Nresids.FL.poly) #Signif; F=3.635, R2=0.25
anova(Cleu.d15Nresids.FL.lm, Cleu.d15Nresids.FL.poly) #poly is better





##BT

SIA %>% filter(Species == "BT") %>% select(d13C) %>% as.matrix() %>% shapiro.test() #non-normal
#see how the resids look after running linear model
SIA %>% filter(Species == "BT") %>% select(d15N) %>% as.matrix() %>% shapiro.test() #normal

SIA %>% filter(Species == "BT") %>% summary() # F=15/M=12; don't need GLS Anova


#test for sex differences

BT.d13C_Sex.lm<- lm(d13C ~ FL*Sex, data = SIA[SIA$Species=="BT",])
plot(BT.d13C_Sex.lm) #meets assumptions
summary(BT.d13C_Sex.lm) #interaction NS; remove

BT.d13C_Sex.lm2<- lm(d13C ~ FL+Sex, data = SIA[SIA$Species=="BT",])
plot(BT.d13C_Sex.lm2) #meets assumptions
summary(BT.d13C_Sex.lm2) #Sex NS; remove

BT.d13C_Sex.lm3<- lm(d13C ~ FL, data = SIA[SIA$Species=="BT",])
plot(BT.d13C_Sex.lm3) #meets assumptions
summary(BT.d13C_Sex.lm3) #FL signif




BT.d15N_Sex.lm<- lm(d15N ~ FL*Sex, data = SIA[SIA$Species=="BT",])
plot(BT.d15N_Sex.lm) #meets assumptions
summary(BT.d15N_Sex.lm) #interaction NS; remove

BT.d15N_Sex.lm2<- lm(d15N ~ FL+Sex, data = SIA[SIA$Species=="BT",])
plot(BT.d15N_Sex.lm2) #meets assumptions
summary(BT.d15N_Sex.lm2) #Sex NS; remove

BT.d15N_Sex.lm3<- lm(d15N ~ FL, data = SIA[SIA$Species=="BT",])
plot(BT.d15N_Sex.lm3) #meets assumptions
summary(BT.d15N_Sex.lm3) #FL signif





#d13C; linear regression fits better, with greater variance at size extremes
ggplot(SIA[SIA$Species=="BT",], aes(x=FL, y=d13C)) + geom_point(size = 3, color = "indianred") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

BT.d13C.FL.lm<- lm(d13C ~ FL, data = SIA[SIA$Species=="BT",])
BT.d13C.FL.poly<- lm(d13C ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="BT",])
summary(BT.d13C.FL.lm) # F=18.35, R2=0.42; Signif
summary(BT.d13C.FL.poly) # F=9.06, R2=0.43; Signif; stick with lm
anova(BT.d13C.FL.lm, BT.d13C.FL.poly)

#check for normality
plot(BT.d13C.FL.lm); shapiro.test(residuals(BT.d13C.FL.lm)) #looks good


BT.d13C.FL.lmResids<- abs(residuals(BT.d13C.FL.lm)) #want to analyze magnitude (abs val) of resids

ggplot(data = SIA[SIA$Species=="BT",], aes(x=FL, y=BT.d13C.FL.lmResids)) + geom_point(size = 3, color = "indianred") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

BT.d13Cresids.FL.lm<- lm(BT.d13C.FL.lmResids ~ FL, data = SIA[SIA$Species=="BT",])
BT.d13Cresids.FL.poly<- lm(BT.d13C.FL.lmResids ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="BT",])
summary(BT.d13Cresids.FL.lm) #NS; F=0.015, R2=0.001
summary(BT.d13Cresids.FL.poly) #Signif; F=6.486, R2=0.35
anova(BT.d13Cresids.FL.lm, BT.d13Cresids.FL.poly) #poly is better


#d15N; quadratic regression fits better, but no change over ontogeny
ggplot(SIA[SIA$Species=="BT",], aes(x=FL, y=d15N)) + geom_point(size = 3, color = "indianred") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

BT.d15N.FL.lm<- lm(d15N ~ FL, data = SIA[SIA$Species=="BT",])
BT.d15N.FL.poly<- lm(d15N ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="BT",])
summary(BT.d15N.FL.lm) # F=8.105, R2=0.24; Signif
summary(BT.d15N.FL.poly) # F=6.842, R2=0.36; lower resid SE, so use this model; Signif
anova(BT.d15N.FL.lm, BT.d15N.FL.poly) #poly is better


BT.d15N.FL.polyResids<- abs(residuals(BT.d15N.FL.poly)) #want to analyze magnitude (abs val) of resids

ggplot(data = SIA[SIA$Species=="BT",], aes(x=FL, y=BT.d15N.FL.polyResids)) + geom_point(size = 3, color = "indianred") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

BT.d15Nresids.FL.lm<- lm(BT.d15N.FL.polyResids ~ FL, data = SIA[SIA$Species=="BT",])
BT.d15Nresids.FL.poly<- lm(BT.d15N.FL.polyResids ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="BT",])
summary(BT.d15Nresids.FL.lm) #NS; F=0.355, R2=0.01
summary(BT.d15Nresids.FL.poly) #NS; F=0.425, R2=0.03
anova(BT.d15Nresids.FL.lm, BT.d15Nresids.FL.poly) #lm better








##BH

SIA %>% filter(Species == "BH") %>% select(d13C) %>% as.matrix() %>% shapiro.test() #non-normal
#see how the resids look after running linear model
SIA %>% filter(Species == "BH") %>% select(d15N) %>% as.matrix() %>% shapiro.test() #normal

SIA %>% filter(Species == "BH") %>% summary() #F=26/M=9; need to use GLS Anova


#check for sex differences in d13C/d15N over FL

BH.d13C_Sex.gls<- gls(d13C ~ FL*Sex, data = SIA[SIA$Species=="BH",], weights = varIdent(form = ~ 1 | Sex), na.action = "na.omit")
plot(BH.d13C_Sex.gls) #might want to change from lm
qqnorm(BH.d13C_Sex.gls, abline = c(0,1)) #normal
anova(BH.d13C_Sex.gls, type = "marginal") #interation NS; remove

BH.d13C_Sex.gls2<- gls(d13C ~ FL + Sex, data = SIA[SIA$Species=="BH",], weights = varIdent(form = ~ 1 | Sex), na.action = "na.omit")
plot(BH.d13C_Sex.gls2) #might want to change from lm
qqnorm(BH.d13C_Sex.gls2, abline = c(0,1)) #normal
anova(BH.d13C_Sex.gls2, type = "marginal") #Sex NS; remove

BH.d13C_Sex.gls3<- gls(d13C ~ FL, data = SIA[SIA$Species=="BH",], weights = varIdent(form = ~ 1 | Sex), na.action = "na.omit")
plot(BH.d13C_Sex.gls3) #might want to change from lm
qqnorm(BH.d13C_Sex.gls3, abline = c(0,1)) #normal
anova(BH.d13C_Sex.gls3, type = "marginal") #FL signif





BH.d15N_Sex.gls<- gls(d15N ~ FL*Sex, data = SIA[SIA$Species=="BH",], weights = varIdent(form = ~ 1 | Sex), na.action = "na.omit")
plot(BH.d15N_Sex.gls) #HOV
qqnorm(BH.d15N_Sex.gls, abline = c(0,1)) #normal
anova(BH.d15N_Sex.gls, type = "marginal") #interation NS; remove

BH.d15N_Sex.gls2<- gls(d15N ~ FL + Sex, data = SIA[SIA$Species=="BH",], weights = varIdent(form = ~ 1 | Sex), na.action = "na.omit")
plot(BH.d15N_Sex.gls2) #HOV
qqnorm(BH.d15N_Sex.gls2, abline = c(0,1)) #normal
anova(BH.d15N_Sex.gls2, type = "marginal") #Sex NS; remove

BH.d15N_Sex.gls3<- gls(d15N ~ FL, data = SIA[SIA$Species=="BH",], weights = varIdent(form = ~ 1 | Sex), na.action = "na.omit")
plot(BH.d15N_Sex.gls3) #might want to change from lm
qqnorm(BH.d15N_Sex.gls3, abline = c(0,1)) #normal
anova(BH.d15N_Sex.gls3, type = "marginal") #FL NS w/ lm






#d13C; quadratic regression fits better, variance decreases over ontogeny and reaches asymptote
ggplot(SIA[SIA$Species=="BH",], aes(x=FL, y=d13C)) + geom_point(size = 3, color = "darkseagreen") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

BH.d13C.FL.lm<- lm(d13C ~ FL, data = SIA[SIA$Species=="BH",])
BH.d13C.FL.poly<- lm(d13C ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="BH",])
summary(BH.d13C.FL.lm) # F=13.84, R2=0.30, Signif
summary(BH.d13C.FL.poly) # F=19.09, R2=0.54, Signif
anova(BH.d13C.FL.lm, BH.d13C.FL.poly) #poly is better

#check for normality
plot(BH.d13C.FL.poly); shapiro.test(residuals(BH.d13C.FL.poly)) #looks good


BH.d13C.FL.polyResids<- abs(residuals(BH.d13C.FL.poly)) #want to analyze magnitude (abs val) of resids

ggplot(data = SIA[SIA$Species=="BH",], aes(x=FL, y=BH.d13C.FL.polyResids)) + geom_point(size = 3, color = "darkseagreen") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

BH.d13Cresids.FL.lm<- lm(BH.d13C.FL.polyResids ~ FL, data = SIA[SIA$Species=="BH",])
BH.d13Cresids.FL.poly<- lm(BH.d13C.FL.polyResids ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="BH",])
summary(BH.d13Cresids.FL.lm) #Signif; F=13.42, R2=0.29
summary(BH.d13Cresids.FL.poly) #Signif; F=19.97, R2=0.56
anova(BH.d13Cresids.FL.lm, BH.d13Cresids.FL.poly) #poly is better


#d15N; linear regression fits better, with decrease in variance over ontogeny
ggplot(SIA[SIA$Species=="BH",], aes(x=FL, y=d15N)) + geom_point(size = 3, color = "darkseagreen") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

BH.d15N.FL.lm<- lm(d15N ~ FL, data = SIA[SIA$Species=="BH",])
BH.d15N.FL.poly<- lm(d15N ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="BH",])
summary(BH.d15N.FL.lm) # F=2.677, R2=0.08, NS
summary(BH.d15N.FL.poly) # F=1.319, R2=0.08, NS; stick w lm
anova(BH.d15N.FL.lm, BH.d15N.FL.poly)


BH.d15N.FL.lmResids<- abs(residuals(BH.d15N.FL.lm)) #want to analyze magnitude (abs val) of resids

ggplot(data = SIA[SIA$Species=="BH",], aes(x=FL, y=BH.d15N.FL.lmResids)) + geom_point(size = 3, color = "darkseagreen") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 2, se=F)

BH.d15Nresids.FL.lm<- lm(BH.d15N.FL.lmResids ~ FL, data = SIA[SIA$Species=="BH",])
BH.d15Nresids.FL.poly<- lm(BH.d15N.FL.lmResids ~ poly(FL, 2, raw = T), data = SIA[SIA$Species=="BH",])
summary(BH.d15Nresids.FL.lm) #Signif; F=9.797, R2=0.23; better model
summary(BH.d15Nresids.FL.poly) #Signif; F=5.001, R2=0.24
anova(BH.d15Nresids.FL.lm, BH.d15Nresids.FL.poly)




### Plot all resid plots together in composite; only signif diff in resids over FL have reg line

cleu.d13C_FL.residsplot<- ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=FL, y=Cleu.d13C.FL.polyResids)) + geom_point(size = 3, color = "steelblue4") + theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), panel.grid = element_blank()) + labs(x = "", y = expression("Abs. values of "*δ^{13}*"C residuals (‰)")) + ylim(0,4)

cleu.d15N_FL.residsplot<- ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=FL, y=Cleu.d15N.FL.polyResids)) + geom_point(size = 3, color = "steelblue4") + theme_bw() + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 1, se=F) + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), panel.grid = element_blank()) + labs(x = "", y = expression("Abs. values of "*δ^{15}*"N residuals (‰)")) + ylim(0,4)


clim.d13C_FL.residsplot<- ggplot(data = SIA[SIA$Species=="BT",], aes(x=FL, y=BT.d13C.FL.lmResids)) + geom_point(size = 3, color = "firebrick") + theme_bw() + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype =1, se=F) + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), panel.grid = element_blank()) + labs(x = "", y = expression("Abs. values of "*δ^{13}*"C residuals (‰)")) + ylim(0,4)

clim.d15N_FL.residsplot<- ggplot(data = SIA[SIA$Species=="BT",], aes(x=FL, y=BT.d15N.FL.polyResids)) + geom_point(size = 3, color = "firebrick") + theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), panel.grid = element_blank()) + labs(x = "", y = expression("Abs. values of "*δ^{15}*"N residuals (‰)")) + ylim(0,4)


stib.d13C_FL.residsplot<- ggplot(data = SIA[SIA$Species=="BH",], aes(x=FL, y=BH.d13C.FL.polyResids)) + geom_point(size = 3, color = "darkseagreen4") + theme_bw() + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", linetype = 1, se=F) + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), panel.grid = element_blank()) + labs(x = "FL (cm)", y = expression("Abs. values of "*δ^{13}*"C residuals (‰)")) + ylim(0,4)

stib.d15N_FL.residsplot<- ggplot(data = SIA[SIA$Species=="BH",], aes(x=FL, y=BH.d15N.FL.lmResids)) + geom_point(size = 3, color = "darkseagreen4") + theme_bw() + geom_smooth(method = "lm", color = "black", se=F) + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), panel.grid = element_blank()) + labs(x = "FL (cm)", y = expression("Abs. values of "*δ^{15}*"N residuals (‰)")) + ylim(0,4)




#Plot all together

SIA.resids.plot<- plot_grid(cleu.d13C_FL.residsplot, NULL, cleu.d15N_FL.residsplot, clim.d13C_FL.residsplot, NULL, clim.d15N_FL.residsplot, stib.d13C_FL.residsplot, NULL, stib.d15N_FL.residsplot, labels = "", nrow = 3, align = "hv", rel_widths = c(1,0.2,1,1,0.2,1,1,0.2,1))

save_plot("SIA_FL residuals plot2.png", SIA.resids.plot, nrow = 3, base_aspect_ratio = 3)





#### Convert Size Class cutoffs from TL to FL; Plot ABF and Stable Isotopes over Ontogeny ####

##Cleu

plot(FL~TL, data = cleu.bfdata)
summary(lm(FL~TL, data = cleu.bfdata)) #R2 = 0.997

#YoY (90 cm TL upper limit)
90*0.824751-3.302571 # 70.93 cm FL

#Juv (160 cm TL upper limit)
160*0.824751-3.302571 #128.66 cm FL

#Sub-adult (210 cm TL upper limit)
210*0.824751-3.302571 #169.90 cm FL


ggplot(cleu.bfdata, aes(x=FL, y=ABF)) + geom_rect(aes(xmin = -Inf, xmax = 70.93, ymin = -Inf, ymax = Inf), fill = "gray90") + geom_rect(aes(xmin = 70.93, xmax = 128.66, ymin = -Inf, ymax = Inf), fill = "gray80") + geom_rect(aes(xmin = 128.66, xmax = 169.9, ymin = -Inf, ymax = Inf), fill = "gray70") + geom_rect(aes(xmin = 169.9, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "gray60") + geom_point(size = 3)





##BT

plot(FL~TL, data = bt.bfdata)
summary(lm(FL~TL, data = bt.bfdata)) #R2 = 0.984

#YoY (83 cm TL upper limit)
83*0.78663+1.50434 # 66.79 cm FL

#Juv (111.5 cm TL upper limit)
111.5*0.78663+1.50434 #89.21 cm FL

#Sub-adult (140 cm TL upper limit)
140*0.78663+1.50434 #111.63 cm FL


ggplot(bt.bfdata, aes(x=FL, y=ABF)) + geom_rect(aes(xmin = -Inf, xmax = 66.79, ymin = -Inf, ymax = Inf), fill = "gray90") + geom_rect(aes(xmin = 66.79, xmax = 89.21, ymin = -Inf, ymax = Inf), fill = "gray80") + geom_rect(aes(xmin = 89.21, xmax = 111.63, ymin = -Inf, ymax = Inf), fill = "gray70") + geom_rect(aes(xmin = 111.63, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "gray60") + geom_point(size = 3)






##BH

plot(FL~TL, data = bh.bfdata)
summary(lm(FL~TL, data = bh.bfdata)) #R2 = 0.988

#YoY (73 cm TL upper limit)
73*0.82114-2.43565 # 57.51 cm FL

#Juv (88.5 cm TL upper limit)
88.5*0.82114-2.43565 #70.24 cm FL


ggplot(bh.bfdata, aes(x=FL, y=ABF)) + geom_rect(aes(xmin = -Inf, xmax = 57.51, ymin = -Inf, ymax = Inf), fill = "gray90") + geom_rect(aes(xmin = 57.51, xmax = 70.24, ymin = -Inf, ymax = Inf), fill = "gray80") + geom_rect(aes(xmin = 70.24, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "gray60") + geom_point(size = 3)


#For PCL
plot(PCL~FL, data = stib.bfdata)
summary(lm(PCL~FL, data = stib.bfdata)) #R2 = 0.983

# y = -1.95778 + 0.94361x







##plot ABF/d13C/d15N over FL (w/ signif regressions for stable isotopes)

cleu.abf.plot<- ggplot(cleu.bfdata, aes(x = FL, y = ABF)) + geom_rect(aes(xmin = -Inf, xmax = 70.93, ymin = -Inf, ymax = Inf), fill = "gray95") + geom_rect(aes(xmin = 70.93, xmax = 128.66, ymin = -Inf, ymax = Inf), fill = "gray87") + geom_rect(aes(xmin = 128.66, xmax = 169.9, ymin = -Inf, ymax = Inf), fill = "gray79") + geom_rect(aes(xmin = 169.9, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "gray71") + geom_point(size = 3, color = "steelblue4") + theme_bw() + labs(x = "", y = "ABF (N)") + theme(axis.text = element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), plot.margin = unit(c(0.5,0.5,0,0.25), "cm")) + scale_x_continuous(limits = c(50,177), breaks = seq(60,180,30))
cleu.d13C.plot<- ggplot(SIA[SIA$Species == "Cleu",], aes(x = FL, y = d13C)) + geom_point(size = 3, color = "steelblue4") + theme_bw() + labs(x = "", y = expression(δ^{13}*"C (‰)")) + theme(axis.text.y = element_text(size=16), axis.text.x = element_blank(), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.margin = unit(c(0.5,0.5,0,0.25), "cm")) + scale_x_continuous(limits = c(50,177), breaks = seq(60,180,30))
cleu.d15N.plot<- ggplot(SIA[SIA$Species == "Cleu",], aes(x = FL, y = d15N)) + geom_point(size = 3, color = "steelblue4") + theme_bw() + labs(x = "", y = expression(δ^{15}*"N (‰)")) + theme(axis.text.y=element_text(size=16), axis.text.x = element_blank(), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)), plot.margin = unit(c(0.5,0.5,0,0.25), "cm")) + scale_x_continuous(limits = c(50,177), breaks = seq(60,180,30))

plot_grid(cleu.abf.plot, cleu.d13C.plot, cleu.d15N.plot, labels = c("A","D","G"), nrow = 3, align = "v")



bt.abf.plot<- ggplot(clim.bfdata, aes(x = FL, y = ABF)) + geom_rect(aes(xmin = -Inf, xmax = 66.79, ymin = -Inf, ymax = Inf), fill = "gray95") + geom_rect(aes(xmin = 66.79, xmax = 89.21, ymin = -Inf, ymax = Inf), fill = "gray87") + geom_rect(aes(xmin = 89.21, xmax = 111.63, ymin = -Inf, ymax = Inf), fill = "gray79") + geom_rect(aes(xmin = 111.63, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "gray71") + geom_point(size = 3, color = "firebrick") + theme_bw() + labs(x = "FL (cm)", y = "") + theme(axis.text = element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) + scale_x_continuous(limits = c(50,137), breaks = seq(50,140,25))
bt.d13C.plot<- ggplot(SIA[SIA$Species == "BT",], aes(x = FL, y = d13C)) + geom_point(size = 3, color = "firebrick") + theme_bw() + labs(x = "", y = "") + theme(axis.text.y = element_text(size=16), axis.text.x = element_blank(), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) + geom_smooth(method = "lm", se = F, color = "black") + scale_x_continuous(limits = c(50,137), breaks = seq(50,140,25))
bt.d15N.plot<- ggplot(SIA[SIA$Species == "BT",], aes(x = FL, y = d15N)) + geom_point(size = 3, color = "firebrick") + theme_bw() + labs(x = "", y = "") + theme(axis.text.y=element_text(size=16), axis.text.x = element_blank(), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", se=F) + scale_x_continuous(limits = c(50,137), breaks = seq(50,140,25))

plot_grid(bt.abf.plot, bt.d13C.plot, bt.d15N.plot, labels = c("B","E","H"), nrow = 3, align = "v")



bh.abf.plot<- ggplot(stib.bfdata, aes(x = FL, y = ABF)) + geom_rect(aes(xmin = -Inf, xmax = 57.51, ymin = -Inf, ymax = Inf), fill = "gray95") + geom_rect(aes(xmin = 57.51, xmax = 70.24, ymin = -Inf, ymax = Inf), fill = "gray87") + geom_rect(aes(xmin = 70.24, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "gray71") + geom_point(size = 3, color = "darkseagreen4") + theme_bw() + labs(x = "", y = "") + theme(axis.text = element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
bh.d13C.plot<- ggplot(SIA[SIA$Species == "BH",], aes(x = FL, y = d13C)) + geom_point(size = 3, color = "darkseagreen4") + theme_bw() + labs(x = "", y = "") + theme(axis.text.y = element_text(size=16), axis.text.x = element_blank(), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) + geom_smooth(method = "lm", formula = y~poly(x, 2, raw = T), color = "black", se=F)
bh.d15N.plot<- ggplot(SIA[SIA$Species == "BH",], aes(x = FL, y = d15N)) + geom_point(size = 3, color = "darkseagreen4") + theme_bw() + labs(x = "", y = "") + theme(axis.text.y=element_text(size=16), axis.text.x = element_blank(), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))

plot_grid(bh.abf.plot, bh.d13C.plot, bh.d15N.plot, labels = c("C","F","I"), nrow = 3, align = "v")


bf_sia.plot<- plot_grid(cleu.d13C.plot, bt.d13C.plot, bh.d13C.plot, cleu.d15N.plot, bt.d15N.plot, bh.d15N.plot, cleu.abf.plot, bt.abf.plot, bh.abf.plot, labels = "", nrow = 3, ncol = 3, align = "hv")
save_plot("ABF_SIA_2.png", bf_sia.plot, nrow = 3, base_aspect_ratio = 5)









#### How do species partition ecological niche space? ####



#ANOVAs to compare mean d13C/d15N

SIA %>% group_by(Species) %>% summarise(n=n()) #unequal sample sizes; weighted GLS ANOVA

d13C.gls<- gls(d13C ~ Species, data = SIA, weights = varIdent(form = ~1 | Species), na.action = "na.omit")
anova(d13C.gls, type = "marginal") #NS
summary(d13C.gls)
plot(d13C.gls) #check for HOV
qqnorm(d13C.gls, abline = c(0,1)) #check normality


d15N.gls<- gls(d15N ~ Species, data = SIA, weights = varIdent(form = ~1 | Species), na.action = "na.omit")
anova(d15N.gls, type = "marginal") #NS
summary(d15N.gls)
plot(d15N.gls) #check for HOV
qqnorm(d15N.gls, abline = c(0,1)) #check normality



library(SIBER)

## load data ##
siber.df<- SIA %>% select(iso1 = d13C, iso2 = d15N, group = Species) %>% mutate(community = rep(1, nrow(SIA)))
siber.data<- createSiberObject(siber.df)


## plot intial data ##

palette(c("cornflowerblue", "indianred", "darkseagreen")) #create palette for points w plot()

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = pchisq(1,2), lty = 1, lwd = 2) #pchisq(1,2) = 0.3935 and is used by Andrew Jackson for estimating SEA (contains ~ 40% of data)
group.hull.args      <- list(lty = 2, col = "gray20")

plotSiberObject(siber.data,
                ax.pad = 1,
                hulls = F, community.hulls.args,
                ellipses = T, group.ellipses.args, #fits SEA (40% probability)
                group.hulls = T, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'),
                points.order = 16
                )


##ggplot2 alternative

#create df for convex hull
conv_hull<- SIA %>% group_by(Species) %>% slice(chull(d13C, d15N))

#use stat_ellipse to create ellipses from MVN distrib

ggplot(data = SIA, aes(x=d13C, y=d15N, color=Species)) + geom_point(size = 3) + theme_bw() + labs(x = expression({delta}^13*C~'(\u2030)'), y = expression({delta}^15*N~'(\u2030)')) + geom_polygon(data = conv_hull, linetype = 2, fill = NA) + stat_ellipse(data = SIA, aes(x=d13C, y=d15N), type = "norm", level = pchisq(1,2), lwd = 1.25) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), legend.position = c(0.88,0.82), legend.background = element_blank(), legend.title = element_blank(), legend.key.size = unit(1.5, 'lines')) + scale_color_manual(values = c('steelblue4','firebrick','darkseagreen4'), labels = c("Bull", "Blacktip", "Bonnethead")) + guides(color = guide_legend(override.aes = list(shape=16, linetype=0, size=4), label.theme = element_text(size = 10)))

ggsave("TA_SEA isotopic niche plot.png", width = 6, height = 4, units = "in", dpi = 600)



##Make plot of only C. leucas for demonstration of niche metrics

#for CR, NR, CD, TA, SEAc, niche var.
ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=d13C, y=d15N)) + geom_point(size = 3) + theme_bw() + labs(x = expression({delta}^13*C~'(\u2030)'), y = expression({delta}^15*N~'(\u2030)')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + geom_segment(aes(x=min(d13C), y=(max(d15N)+min(d15N))/2, xend=max(d13C), yend=(max(d15N)+min(d15N))/2), color = "blue", size = 1.5, arrow = arrow(length = unit(0.05, "npc"), ends = "both"))

ggsave("CR example plot.png", width = 6, height = 4, units = "in", dpi = 600)



ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=d13C, y=d15N)) + geom_point(size = 3) + theme_bw() + labs(x = expression({delta}^13*C~'(\u2030)'), y = expression({delta}^15*N~'(\u2030)')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + geom_segment(aes(x=(max(d13C)+min(d13C))/2, y=min(d15N), xend=(max(d13C)+min(d13C))/2, yend=max(d15N)), color = "blue", size = 1.5, arrow = arrow(length = unit(0.05, "npc"), ends = "both"))

ggsave("NR example plot.png", width = 6, height = 4, units = "in", dpi = 600)



ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=d13C, y=d15N)) + geom_point(size = 3) + theme_bw() + labs(x = expression({delta}^13*C~'(\u2030)'), y = expression({delta}^15*N~'(\u2030)')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + geom_segment(aes(x=d13C, y=d15N, xend=mean(d13C), yend=mean(d15N)), color = "blue", size = 0.5, arrow = arrow(length = unit(0.03, "npc"), ends = "both")) + geom_point(aes(x=mean(d13C), y=mean(d15N)), color = "red", size=4)

ggsave("CD example plot.png", width = 6, height = 4, units = "in", dpi = 600)



conv_hull_cleu<- SIA[SIA$Species=="Cleu",] %>% slice(chull(d13C, d15N))

ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=d13C, y=d15N)) + geom_point(size = 3) + theme_bw() + labs(x = expression({delta}^13*C~'(\u2030)'), y = expression({delta}^15*N~'(\u2030)')) + geom_polygon(data = conv_hull_cleu, size = 1, linetype = 1, fill = NA, color="blue") + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

ggsave("TA example plot.png", width = 6, height = 4, units = "in", dpi = 600)



ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=d13C, y=d15N)) + geom_point(size = 3) + theme_bw() + labs(x = expression({delta}^13*C~'(\u2030)'), y = expression({delta}^15*N~'(\u2030)')) + stat_ellipse(type = "norm", level = pchisq(1,2), lwd = 1.25, color = "blue") + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

ggsave("SEAc example plot.png", width = 6, height = 4, units = "in", dpi = 600)


BH.SIA.dat<- SIA[SIA$Species=="BH",]
BH.SIA.dat$predicted<- predict(BH.d15N.FL.lm)

ggplot(data = BH.SIA.dat, aes(x=FL, y=d15N)) + geom_point(size = 3) + theme_bw() + labs(x = "FL (cm)", y = expression({delta}^15*N~'(\u2030)')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + stat_smooth(method = "lm", se=F) + geom_segment(aes(xend = FL, yend = predicted), arrow = arrow(length = unit(0.03, "npc"), ends = "both"), color = "red")

ggsave("Niche var plot.png", width = 6, height = 4, units = "in", dpi = 600)



ggplot(data = BH.SIA.dat, aes(x=FL, y=BH.d15N.FL.lmResids)) + geom_point(size = 3) + theme_bw() + geom_smooth(method = "lm", se=F) + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), panel.grid = element_blank()) + labs(x = "FL (cm)", y = expression("Abs. values of "*δ^{15}*"N residuals (‰)"))

ggsave("Niche var plot2.png", width = 6, height = 4, units = "in", dpi = 600)




## explore niche width (TA, SEA, SEAc) ##
groupMetricsML(siber.data) #Cleu > BH > BT

## explore Layman's metrics ##

laymanMetrics(SIA[SIA$Species=="Cleu","d13C"], SIA[SIA$Species=="Cleu","d15N"])$metrics #Cleu
laymanMetrics(SIA[SIA$Species=="BT","d13C"], SIA[SIA$Species=="BT","d15N"])$metrics #BT
laymanMetrics(SIA[SIA$Species=="BH","d13C"], SIA[SIA$Species=="BH","d15N"])$metrics #BH




## run Bayesian model ##

# options for running jags
params <- list()
params$n.iter <- 2 * 10^4   # number of iterations to run the model for
params$n.burnin <- 1 * 10^3 # discard the first set of values
params$n.thin <- 10     # thin the posterior by this many
params$n.chains <- 2        # run this many chains

# set save.output = TRUE
params$save.output = TRUE
params$save.dir = tempdir()

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior (conjugate prior of MVN)
# on the covariance matrix Sigma, and a vague normal prior on the
# means. Fitting is via the JAGS method.

#First, test each species for multivariate normality (Royston's H-test via Shapiro-Wilk's W)
library(MVN)
SIA %>% filter(Species == "Cleu") %>% dplyr::select(d13C, d15N) %>% mvn(mvnTest = "royston") #good
SIA %>% filter(Species == "BT") %>% dplyr::select(d13C, d15N) %>% mvn(mvnTest = "royston") #d13C not normal
SIA %>% filter(Species == "BH") %>% dplyr::select(d13C, d15N) %>% mvn(mvnTest = "royston") #d13C not normal


ellipses.posterior <- siberMVN(siber.data, params, priors)



## test convergence ##
library(coda)

# get a list of all the files in the save directory
all.files <- dir(params$save.dir, full.names = TRUE)

# find which ones are jags model files
model.files <- all.files[grep("jags_output", all.files)]



# test convergence for group one
group1 <- 1

load(model.files[group1])

plot(output) #check trace plots
gelman.diag(output, multivariate = FALSE)
gelman.plot(output, auto.layout = FALSE)


# test convergence for group two
group2 <- 2

load(model.files[group2])

plot(output)
gelman.diag(output, multivariate = FALSE)
gelman.plot(output, auto.layout = FALSE)


# test convergence for group two
group3 <- 3

load(model.files[group3])

plot(output)
gelman.diag(output, multivariate = FALSE)
gelman.plot(output, auto.layout = FALSE)




# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

png("SEA_B niche size plot.png", width = 6, height = 4, units = "in", res = 600)

par(mar = c(4,6,3,2) + 0.1) #adjust margins to fit y-axis title
siberDensityPlot(SEA.B, xticklabels = c("","",""),
                 xlab = "",
                 ylab = "",
                 las = 1,
                 ct = "mode", #black dot represents mode
                 prn = T,
                 ylab.line = 3.5,
                 ylims = c(0,10),
                 axes = F
                 )
title(ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ), cex.lab = 1.5)
axis(1, at = 1:3, labels = c("Bull","Blacktip","Bonnethead"), cex.axis = 1.5)
axis(2, cex.axis = 1.25)
box() #manually add box around plot

dev.off()

par(mar = c(5,4,4,2) + 0.1) #reset back to original margins


#alternative densityplot version
library(reshape2)
SEA.B.melt<- melt(SEA.B, varnames = names(SEA.B))
SEA.B.melt$Var2<- factor(SEA.B.melt$Var2)
levels(SEA.B.melt$Var2)<- c("Cleu", "BT", "BH")

ggplot(SEA.B.melt, aes(fill = Var2, x = value)) + geom_density(alpha = .5)




## Test probability of one SEA.B being larger than another (BEST package if needed)

#Prob that Cleu SEA.B is larger than BT SEA.B (100%)
NROW(SEA.B[SEA.B[,1] > max(SEA.B[,2]),1]) / 4000 * 100

#Prob that BH SEA.B is larger than BT SEA.B (99.63%)
NROW(SEA.B[SEA.B[,3] > max(SEA.B[,2]),3]) / 4000 * 100

#Prob that Cleu SEA.B is larger than BH SEA.B (99.48%)
NROW(SEA.B[SEA.B[,1] > max(SEA.B[,3]),1]) / 4000 * 100



## measure niche overlap ##

#TA

cleu.hull<- siberConvexhull(SIA[SIA$Species=="Cleu",9], SIA[SIA$Species=="Cleu",8])
bt.hull<- siberConvexhull(SIA[SIA$Species=="BT",9], SIA[SIA$Species=="BT",8])
bh.hull<- siberConvexhull(SIA[SIA$Species=="BH",9], SIA[SIA$Species=="BH",8])

library(spatstat.utils) #for overlap.xypolygon fxn

Cleu_BT.overlap.TA<- overlap.xypolygon(list(x = cleu.hull$xcoords, y = cleu.hull$ycoords), list(x = bt.hull$xcoords, y = bt.hull$ycoords)) #as area (per mil ^ 2)
#measure percentage of overlap directionally (i.e. Cleu onto BT and BT onto Cleu)
Cleu.onto.BT.TA<- Cleu_BT.overlap.TA / as.numeric(bt.hull$TA) * 100
BT.onto.Cleu.TA<- Cleu_BT.overlap.TA / as.numeric(cleu.hull$TA) * 100


BT_BH.overlap.TA<- overlap.xypolygon(list(x = bt.hull$xcoords, y = bt.hull$ycoords), list(x = bh.hull$xcoords, y = bh.hull$ycoords)) #as area (per mil ^ 2)
#measure percentage of overlap directionally (i.e. BT onto BH and BH onto BT)
BT.onto.BH.TA<- BT_BH.overlap.TA / as.numeric(bh.hull$TA) * 100
BH.onto.BT.TA<- BT_BH.overlap.TA / as.numeric(bt.hull$TA) * 100


Cleu_BH.overlap.TA<- overlap.xypolygon(list(x = cleu.hull$xcoords, y = cleu.hull$ycoords), list(x = bh.hull$xcoords, y = bh.hull$ycoords)) #as area (per mil ^ 2)
#measure percentage of overlap directionally (i.e. Cleu onto BH and BH onto Cleu)
Cleu.onto.BH.TA<- Cleu_BH.overlap.TA / as.numeric(bh.hull$TA) * 100
BH.onto.Cleu.TA<- Cleu_BH.overlap.TA / as.numeric(cleu.hull$TA) * 100



#SEAc

Cleu_BT.overlap.SEAc<- maxLikOverlap("1.Cleu", "1.BT", siber.data, p.interval = NULL, n = 100) #for ~40% prediction interval of SEAc
#measure percentage of overlap directionally (i.e. Cleu onto BT and BT onto Cleu)
Cleu.onto.BT.SEAc<- Cleu_BT.overlap.SEAc[3] / Cleu_BT.overlap.SEAc[2] * 100
BT.onto.Cleu.SEAc<- Cleu_BT.overlap.SEAc[3] / Cleu_BT.overlap.SEAc[1] * 100


BT_BH.overlap.SEAc<- maxLikOverlap("1.BT", "1.BH", siber.data, p.interval = NULL, n = 100)
#measure percentage of overlap directionally (i.e. BT onto BH and BH onto BT)
BT.onto.BH.SEAc<- BT_BH.overlap.SEAc[3] / BT_BH.overlap.SEAc[2] * 100
BH.onto.BT.SEAc<- BT_BH.overlap.SEAc[3] / BT_BH.overlap.SEAc[1] * 100


Cleu_BH.overlap.SEAc<- maxLikOverlap("1.Cleu", "1.BH", siber.data, p.interval = NULL, n = 100)
#measure percentage of overlap directionally (i.e. Cleu onto BH and BH onto Cleu)
Cleu.onto.BH.SEAc<- Cleu_BH.overlap.SEAc[3] / Cleu_BH.overlap.SEAc[2] * 100
BH.onto.Cleu.SEAc<- Cleu_BH.overlap.SEAc[3] / Cleu_BH.overlap.SEAc[1] * 100



#SEAb via nicheROVER for CIs

library(nicheROVER)
library(viridis)

aggregate(SIA[c("d13C","d15N")], SIA[2], mean) #view SI means by Species

set.seed(123)

nsamples <- 1000
SEAb.par <- tapply(1:nrow(SIA), SIA$Species, function(ii) niw.post(nsamples = nsamples, X = SIA[ii, c("d13C","d15N")]))

niche.par.plot(SEAb.par, col = c(viridis(3)), plot.mu = T, plot.Sigma = F) #view distrib of means for d13C and d15N
legend("topleft", legend = names(SEAb.par), fill = viridis(3), bty = "n")

nsamples <- 10
SEAb.par2 <- tapply(1:nrow(SIA), SIA$Species, function(ii) niw.post(nsamples = nsamples, X = SIA[ii, c("d13C","d15N")]))
# format data for plotting function
SIA.data <- tapply(1:nrow(SIA), SIA$Species, function(ii) X = SIA[ii, c("d13C","d15N")])

niche.plot(niche.par = SEAb.par2, niche.data = SIA.data, alpha = 0.40, pfrac = 0.1, iso.names = expression(delta^{
  15
} * N, delta^{
  13
} * C), col = viridis(3), xlab = expression("Isotope Ratio (per mil)")) #using 40% of niche region


# niche overlap plots for 40% niche region sizes

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher
# accuracy.  the variable over.stat can be supplied directly to the
# overlap.plot function

over.stat.40 <- overlap(SEAb.par, nreps = 10000, nprob = 10000, alpha = pchisq(1,2))

# The mean overlap metrics calculated across iteratations for the niche
# region size (alpha = 0.40) can be calculated and displayed
# in an array.
over.mean.40 <- apply(over.stat.40, c(1:2), mean) * 100
round(over.mean.40, 2)

over.cred.40 <- apply(over.stat.40 * 100, c(1:2), quantile, prob = c(0.025, 0.975), na.rm = TRUE)
round(over.cred.40, 2)  # display 95% credible intervals for alpha = 0.40 niche region


overlap.plot(over.stat.40, col = viridis(3), mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 40%")







# niche overlap plots for 95% niche region sizes

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher
# accuracy.  the variable over.stat can be supplied directly to the
# overlap.plot function

over.stat.95 <- overlap(SEAb.par, nreps = 10000, nprob = 10000, alpha = 0.95)

# The mean overlap metrics calculated across iteratations for the niche
# region size (alpha = 0.95) can be calculated and displayed
# in an array.
over.mean.95 <- apply(over.stat.95, c(1:2), mean) * 100
round(over.mean.95, 2)

over.cred.95 <- apply(over.stat.95 * 100, c(1:2), quantile, prob = c(0.025, 0.975), na.rm = TRUE)
round(over.cred.95, 2)  # display 95% credible intervals for alpha = 0.95 niche region


overlap.plot(over.stat.95, col = viridis(3), mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")








#overlap expressed as area (0/00^2) via SIBER

Cleu_BT.overlap.area<- bayesianOverlap(ellipse1 = "1.BT", ellipse2 = "1.Cleu", ellipses.posterior = ellipses.posterior, draws = 1000, p.interval = 0.95, do.plot = T)
#overlapped expressed as percentage of non-overlapping area
Cleu_BT.overlap.per<- (Cleu_BT.overlap.area[,3] / (Cleu_BT.overlap.area[,2] + Cleu_BT.overlap.area[,1] - Cleu_BT.overlap.area[,3])) * 100


BH_BT.overlap.area<- bayesianOverlap(ellipse1 = "1.BT", ellipse2 = "1.BH", ellipses.posterior = ellipses.posterior, draws = 100, p.interval = 0.95, do.plot = T)
BH_BT.overlap.per<- (BH_BT.overlap.area[,3] / (BH_BT.overlap.area[,2] + BH_BT.overlap.area[,1] - BH_BT.overlap.area[,3])) * 100


Cleu_BH.overlap.area<- bayesianOverlap(ellipse1 = "1.BH", ellipse2 = "1.Cleu", ellipses.posterior = ellipses.posterior, draws = 100, p.interval = 0.95, do.plot = T)
Cleu_BH.overlap.per<- (Cleu_BH.overlap.area[,3] / (Cleu_BH.overlap.area[,2] + Cleu_BH.overlap.area[,1] - Cleu_BH.overlap.area[,3])) * 100





###=============================================================================================###

#### Plot Bayesian posterior ellipses as example for C. leucas ####

# how many of the posterior draws do you want?
n.posts <- 20

# for a standard ellipse use
p.ell <- pchisq(1,2)


# a list to store the results
all_ellipses <- list()

# loop over groups
for (i in 1:length(ellipses.posterior)){

  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL

  for ( j in 1:n.posts){

    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)

    # mean
    mu     <- ellipses.posterior[[i]][j,5:6]

    # ellipse points

    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)


    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))

  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}

ellipse_df <- data.frame(all_ellipses[1])



ellipse_df <- dplyr::rename(ellipse_df, d13C = x, d15N = y)


ggplot(data = SIA[SIA$Species=="Cleu",], aes(x=d13C, y=d15N)) + geom_point(size = 3) + theme_bw() + labs(x = expression({delta}^13*C~'(\u2030)'), y = expression({delta}^15*N~'(\u2030)')) + geom_path(data = ellipse_df, aes(group = rep), alpha = 0.3, color = "blue") + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

ggsave("SEAb example plot.png", width = 6, height = 4, units = "in", dpi = 600)






#### Directly correlate f'(x) of RMA and Bayesian regressions w/ isotope values over ontogeny ####


### Cleu

## BF funs

cleu.small.RMA<- function(x) (10^-3.524023)*(x^2.997809)
cleu.small.Bayes<- function(x) (0.0001431931)*(x^3.16)

cleu.large.RMA<- function(x) (10^-0.77765161)*(x^1.680699)
cleu.large.Bayes<- function(x) (0.2630845)*(x^1.59)


## Find slope (f'(x)) of funs

cleu.small.RMA_deriv<- Deriv(cleu.small.RMA, "x", nderiv = 1)
cleu.small.Bayes_deriv<- Deriv(cleu.small.Bayes, "x", nderiv = 1)

cleu.large.RMA_deriv<- Deriv(cleu.large.RMA, "x", nderiv = 1)
cleu.large.Bayes_deriv<- Deriv(cleu.large.Bayes, "x", nderiv = 1)




### Clim

## BF funs

clim.small.RMA<- function(x) (10^-3.617967)*(x^2.895455)
clim.small.Bayes<- function(x) (2.896163e-06)*(x^3.87)

clim.large.RMA<- function(x) (10^-3.066041)*(x^2.635060)
clim.large.Bayes<- function(x) (0.02280475)*(x^1.95)


## Find slope (f'(x)) of funs

clim.small.RMA_deriv<- Deriv(clim.small.RMA, "x", nderiv = 1)
clim.small.Bayes_deriv<- Deriv(clim.small.Bayes, "x", nderiv = 1)

clim.large.RMA_deriv<- Deriv(clim.large.RMA, "x", nderiv = 1)
clim.large.Bayes_deriv<- Deriv(clim.large.Bayes, "x", nderiv = 1)




### Stib

## BF funs

stib.small.RMA<- function(x) (10^-4.484093)*(x^3.198358)
stib.small.Bayes<- function(x) (2.730506e-07)*(x^4.32)

stib.large.RMA<- function(x) (10^-1.277219)*(x^1.6110855)
stib.large.Bayes<- function(x) (2.903973)*(x^0.72)


## Find slope (f'(x)) of funs

stib.small.RMA_deriv<- Deriv(stib.small.RMA, "x", nderiv = 1)
stib.small.Bayes_deriv<- Deriv(stib.small.Bayes, "x", nderiv = 1)

stib.large.RMA_deriv<- Deriv(stib.large.RMA, "x", nderiv = 1)
stib.large.Bayes_deriv<- Deriv(stib.large.Bayes, "x", nderiv = 1)



####============================================================================================####

#### Initial Code for Compiling SIA data ####


t1<- read.csv("Tray 1_7.5.16.csv", header = T, sep = ",")
t2<- read.csv("Tray 2_7.18.16.csv", header = T, sep = ",")
t3<- read.csv("Tray 1_2.20.18.csv", header = T, sep = ",")

SIA<- rbind(t1,t2,t3)
names(SIA)[1]<- "SharkID"

SIA$SharkID<- gsub("/","-",SIA$SharkID)

#Check C:N
SIA$C.N<- SIA$wt.per.C/SIA$wt.per.N #calc C:N for each sample
#do any samples have a C:N > 3.5?
SIA$C.N[SIA$C.N >= 3.5] #no, all samples < 3.5




##Import and join TL/BF data##
setwd("~/Documents/Shark Project/BF Data")

cleu<- read.csv("Cleu BF Values.csv", header = T, sep = ",")
bt<- read.csv("BT BF Values.csv", header = T, sep = ",")
bh<- read.csv("BH BF Values.csv", header = T, sep = ",")

cleu<- cleu[,1:14]
bt<- bt[,1:14]
bh<- bh[,1:14]
bfdata<- rbind(cleu,bt,bh) #bring all data into single df

SIA2<- left_join(SIA,bfdata,by="SharkID")


#Subset by species
bh.locs<- grep(pattern='BH',x=SIA2$SharkID)
bh.sia<- SIA2[bh.locs,]
bh.sia$Species<- rep("BH", length(bh.sia$SharkID))

bt.locs<- grep(pattern = "BT", x=SIA2$SharkID)
bt.sia<- SIA2[bt.locs,]
bt.sia$Species<- rep("BT", length(bt.sia$SharkID))

cleu.locs<- grep(pattern = "Cleu", x= SIA2$SharkID)
cleu.sia<- SIA2[cleu.locs,]
cleu.sia$Species<- rep("Cleu", length(cleu.sia$SharkID))

SIA3<- rbind(cleu.sia,bt.sia,bh.sia)
SIA3<- SIA3 %>% select(SharkID, Species, TL, Age, everything())


#write master .csv file
setwd("~/Documents/Shark Project/Diet/SIA")
write.csv(SIA3, "SIA Master.csv", row.names = F)










##### Code for manipulating Tox SIA data after import #####

SIA.tox<- read.csv("Tray 1_2.20.18.csv", header = T, sep = ",")

names(SIA.tox)[1]<- "SharkID"


#Check C:N
avg.per.N<- mean(SIA.tox$wt.per.N) #find avg %N
avg.per.C<- mean(SIA.tox$wt.per.C) #find avg %C
C.N<- avg.per.C/avg.per.N




#Subset by species
bh.locs<- grep(pattern='BH',x=SIA.tox$SharkID)
bh.sia<- SIA.tox[bh.locs,]
bh.sia$species<- rep("BH", length(bh.sia$SharkID))

bt.locs<- grep(pattern = "BT", x=SIA.tox$SharkID)
bt.sia<- SIA.tox[bt.locs,]
bt.sia$species<- rep("BT", length(bt.sia$SharkID))

cleu.locs<- grep(pattern = "Cleu", x= SIA.tox$SharkID)
cleu.sia<- SIA.tox[cleu.locs,]
cleu.sia$species<- rep("Cleu", length(cleu.sia$SharkID))

SIA.tox2<- rbind(cleu.sia,bt.sia,bh.sia)
names(liverWW)[1]<- "SharkID"
SIA.tox3<- left_join(SIA.tox2,liverWW, by = "SharkID")

SIA.tox3[1,8]<- "Cleu"
SIA.tox3[1,11]<- 196.3

#Make biplots and look at trends over TL
ggplot(SIA.tox3[SIA.tox3$species=="BH",], aes(d13C,d15N)) + geom_point(size=3)
ggplot(SIA.tox3[SIA.tox3$species=="BH",], aes(TL,d15N)) + geom_point(size=3)

ggplot(SIA.tox3[SIA.tox3$species=="BT",], aes(d13C,d15N)) + geom_point(size=3)
ggplot(SIA.tox3[SIA.tox3$species=="BT",], aes(TL,d15N)) + geom_point(size=3)

ggplot(SIA.tox3[SIA.tox3$species=="Cleu",], aes(d13C,d15N)) + geom_point(size=3)
ggplot(SIA.tox3[SIA.tox3$species=="Cleu",], aes(TL,d15N)) + geom_point(size=3)


#Look at all 3 species together
ggplot(SIA.tox3, aes(d13C,d15N,color=species)) + geom_point(size=3) + theme_bw()
ggplot(SIA.tox3, aes(TL,d15N,color=species)) + geom_point(size=3) + geom_smooth(method = "lm", se=T) + theme_bw()

#overlapped w/ old data
ggplot(SIA.tox3, aes(d13C,d15N,color=species)) + geom_point(size=4, alpha = 0.75) + geom_point(data = SIA3, aes(d13C,d15N))

SIA3.a<- dplyr::select(SIA3, c(species,TL,d13C,d15N))
SIA.tox3a<- dplyr::select(SIA.tox3, c(species,TL,d13C,d15N))
ggplot(rbind(SIA3.a,SIA.tox3a), aes(TL,d15N,color=species)) + geom_point(size=3) + geom_smooth(method = "lm", se=T)
