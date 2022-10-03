###Importing all .csv files for data manipulation###

setwd("~/Documents/Shark Project/Digitized Points/")
temp = list.files(pattern="*.csv")
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
         read.csv), envir = .GlobalEnv)



sums<- function(x) {colSums(select(x, A.CSA, Mass), na.rm=T)}

##Couldn't reshape the ACSA and Mass data in R; entered in BF csv files##

library(MASS) #for exploratory stepAIC(x, direction='both')
library( MuMIn) #for further analyses
library(rpart) #for using regression trees in rpart()
library(randomForest) #for using randomForest w perms

setwd("~/Documents/Shark Project/BF Data")
bh<- read.csv('BH BF Values.csv', header = T, sep = ',')
bt<- read.csv('BT BF Values.csv', header = T, sep = ',')
cleu<- read.csv('Cleu BF Values.csv', header = T, sep = ',')

acsa<- na.omit(cleu$ACSA) #since Cleu_Ara16 is missing raw files
mass<- na.omit(cleu$Mass)
ama<- cleu$AMA[c(1:14, 16:20)]
abf<- cleu$ABF[c(1:14,16:20)]


cor(bh$AMA, bh$ACSA) #determine if vars are independent
cor(bh$AMA, bh$Mass)
cor(bh$ACSA, bh$Mass)
cor(bt$AMA, bt$ACSA)
cor(bt$AMA, bt$Mass)
cor(bt$ACSA, bt$Mass)
cor(acsa, mass)
cor(acsa, ama)
cor(mass, ama)
#there appears to be high corr for all ACSA-Mass relationships

#stepwise regression using AIC
bh.glm<- glm(bh$ABF ~ bh$AMA + bh$ACSA + bh$Mass)
bh.step<- stepAIC(bh.glm, direction = 'both')
bh.step$anova
bt.glm<- glm(bt$ABF ~ bt$AMA + bt$ACSA + bt$Mass)
bt.step<- stepAIC(bt.glm, direction = 'both')
bt.step$anova
cleu.glm<- glm(abf ~ ama + acsa + mass)
cleu.step<- stepAIC(cleu.glm, direction = 'both')
cleu.step$anova

#adaptive regression by mixing arm.glm; uses AIC
bh.arm<- arm.glm(bh.glm, R=500, weight.by = 'aic')
summary(bh.arm)
bt.arm<- arm.glm(bt.glm, R=500, weight.by = 'aic')
summary(bt.arm)
cleu.arm<- arm.glm(cleu.glm, R=500, weight.by = 'aic')
summary(cleu.arm)

#comparison of models of AICc values using model.sel()
bh.glm2<- update(bh.glm, .~. -bh$Mass)
bh.glm3<- update(bh.glm, .~. -bh$ACSA)
bh.glm4<- update(bh.glm, .~. -bh$AMA)
bhAICc<- model.sel(bh.glm, bh.glm2, bh.glm3, bh.glm4)
importance(bhAICc)

bt.glm2<- update(bt.glm, .~. -bt$Mass)
bt.glm3<- update(bt.glm, .~. -bt$ACSA)
bt.glm4<- update(bt.glm, .~. -bt$AMA)
btAICc<- model.sel(bt.glm, bt.glm2, bt.glm3, bt.glm4)
importance(btAICc)

cleu.glm2<- update(cleu.glm, .~. -mass)
cleu.glm3<- update(cleu.glm, .~. -acsa)
cleu.glm4<- update(cleu.glm, .~. -ama)
cleuAICc<- model.sel(cleu.glm, cleu.glm2, cleu.glm3, cleu.glm4)
importance(cleuAICc)


#Inclusion of interaction term (ACSA*Mass) in model
bh.glm5<- glm(bh$ABF~bh$AMA+bh$ACSA+bh$Mass+bh$ACSA*bh$Mass)
bh.step5<- stepAIC(bh.glm5, direction = 'both')
bh.step5$anova
bhAICc2<- model.sel(bh.glm, bh.glm2, bh.glm3, bh.glm4, bh.glm5)
importance(bhAICc2)

bt.glm5<- glm(bt$ABF~bt$AMA+bt$ACSA+bt$Mass+bt$ACSA*bt$Mass)
bt.step5<- stepAIC(bt.glm5, direction = 'both')
bt.step5$anova
btAICc2<- model.sel(bt.glm, bt.glm2, bt.glm3, bt.glm4, bt.glm5)
importance(btAICc2)

cleu.glm5<- glm(abf~ama+acsa+mass+acsa*mass)
cleu.step5<- stepAIC(cleu.glm5, direction = 'both')
cleu.step5$anova
cleuAICc2<- model.sel(cleu.glm, cleu.glm2, cleu.glm3, cleu.glm4, cleu.glm5)
importance(cleuAICc2)


#Make viewing of model.sel output easier with model.avg output
bh.table<- as.data.frame(bhAICc2)[,6:10]
bh.table[,2:3]<- round(bh.table[,2:3], 2)
bh.table[,4:5]<- round(bh.table[,4:5], 3)
bh.table$Model<- rownames(bh.table)
for(i in 1:nrow(bh.table)) bh.table$Model[i]<- as.character(formula(paste(bh.table$Model[i])))[3]
bh.table<- bh.table[,c(6,1:5)]
summary(model.avg(bhAICc2))

bt.table<- as.data.frame(btAICc2)[,6:10]
bt.table[,2:3]<- round(bt.table[,2:3], 2)
bt.table[,4:5]<- round(bt.table[,4:5], 3)
bt.table$Model<- rownames(bt.table)
for(i in 1:nrow(bt.table)) bt.table$Model[i]<- as.character(formula(paste(bt.table$Model[i])))[3]
bt.table<- bt.table[,c(6,1:5)]
summary(model.avg(btAICc2))

cleu.table<- as.data.frame(cleuAICc2)[,6:10]
cleu.table[,2:3]<- round(cleu.table[,2:3], 2)
cleu.table[,4:5]<- round(cleu.table[,4:5], 3)
cleu.table$Model<- rownames(cleu.table)
for(i in 1:nrow(cleu.table)) cleu.table$Model[i]<- as.character(formula(paste(cleu.table$Model[i])))[3]
cleu.table<- cleu.table[,c(6,1:5)]
summary(model.avg(cleuAICc2))
