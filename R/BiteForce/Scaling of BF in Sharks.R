
### Estimate polynomial regressions, find where 2nd derivative = 0, and compute separate scaling curves for each species ###

library(dplyr)
library(ggplot2)
library(car)
library(Deriv)
library(R2jags)
library(doParallel)
library(broom)
library(DHARMa)
library(AICcmodavg)


setwd("~/Documents/Shark Project/BF Data")


#### Polynomial Regression ####

## C. leucas ##

cleu.bfdata<- read.csv("Cleu BF Values.csv", header=TRUE, sep=",")


#create models
cleu.bf.mod1<- lm(ABF ~ FL, data = cleu.bfdata)
cleu.bf.mod2<- lm(ABF ~ poly(FL, 2, raw = T), data = cleu.bfdata) #quadratic term
cleu.bf.mod3<- lm(ABF ~ poly(FL, 3, raw = T), data = cleu.bfdata) #cubic term
cleu.bf.mod4<- lm(ABF ~ poly(FL, 4, raw = T), data = cleu.bfdata) #quartic term

cleu.bf.mods<- list("mod1" = cleu.bf.mod1, "mod2" = cleu.bf.mod2, "mod3" = cleu.bf.mod3, "mod4" = cleu.bf.mod4)

#test for differences among models by SS (anova fxn) and AICc
anova(cleu.bf.mod1, cleu.bf.mod2, cleu.bf.mod3, cleu.bf.mod4) #model 3 is best
aictab(cand.set = cleu.bf.mods, second.ord = T) #model 3 is best

#Investigate the final model (Model 3)
summary(cleu.bf.mod3) #provides coefficients for model, p-value for overall model, and R2


#Plot polynomial regression from each model for comparison

seq.x = seq(min(cleu.bfdata$FL), max(cleu.bfdata$FL), length.out = 100) #create sequence within FL range
pred.y.mod1 = predict(cleu.bf.mod1, data.frame(FL = seq.x)) #based on FL sequence, estimate ABF for model 1
pred.y.mod2 = predict(cleu.bf.mod2, data.frame(FL = seq.x, FL2 = seq.x^2)) #estimate for model 2
pred.y.mod3 = predict(cleu.bf.mod3, data.frame(FL = seq.x, FL2 = seq.x^2, FL3 = seq.x^3)) #estimate for model 3
pred.y.mod4 = predict(cleu.bf.mod4, data.frame(FL = seq.x, FL2 = seq.x^2, FL3 = seq.x^3, FL4 = seq.x^4)) #estimate for model 4

plot(ABF ~ FL, data = cleu.bfdata, pch = 16)
lines(seq.x, pred.y.mod1, lty=2, lwd=2, col="red")
lines(seq.x, pred.y.mod2, lty=2, lwd=2, col="orange")
lines(seq.x, pred.y.mod3, lty=1, lwd=3, col="blue")
lines(seq.x, pred.y.mod4, lty=2, lwd=2, col="green")



#Extract coeffs from cleu.bf.model3 and create fxn to find where 2nd derivative = 0

cleu.bf.int<- as.numeric(cleu.bf.mod3$coefficients[1]) #intercept
cleu.bf.FL<- as.numeric(cleu.bf.mod3$coefficients[2]) #FL
cleu.bf.FL2<- as.numeric(cleu.bf.mod3$coefficients[3]) #FL^2
cleu.bf.FL3<- as.numeric(cleu.bf.mod3$coefficients[4]) #FL^3

cleu.bf.fxn<- function(x) cleu.bf.FL3*x^3 + cleu.bf.FL2*x^2 + cleu.bf.FL*x + cleu.bf.int

#plot
ggplot(data = cleu.bfdata, aes(x = FL)) + geom_point(aes(y = ABF), size = 3) + stat_function(fun = cleu.bf.fxn, size = 1, color = "blue") + theme_bw()





## C. limbatus ##

clim.bfdata<- read.csv("BT BF Values.csv", header=TRUE, sep=",")


#create models
clim.bf.mod1<- lm(ABF ~ FL, data = clim.bfdata)
clim.bf.mod2<- lm(ABF ~ poly(FL, 2, raw = T), data = clim.bfdata) #quadratic term
clim.bf.mod3<- lm(ABF ~ poly(FL, 3, raw = T), data = clim.bfdata) #cubic term
clim.bf.mod4<- lm(ABF ~ poly(FL, 4, raw = T), data = clim.bfdata) #quartic term

clim.bf.mods<- list("mod1" = clim.bf.mod1, "mod2" = clim.bf.mod2, "mod3" = clim.bf.mod3, "mod4" = clim.bf.mod4)


#test for differences among models by SS (anova fxn) and AICc
anova(clim.bf.mod1, clim.bf.mod2, clim.bf.mod3, clim.bf.mod4) #model 3 is best
aictab(cand.set = clim.bf.mods, second.ord = T) #model 3 dAICc only +1; stay with model 3

#Investigate the final model (Model 3)
summary(clim.bf.mod3) #provides coefficients for model, p-value for overall model, and R2


#Plot polynomial regression from each model for comparison

seq.x = seq(min(clim.bfdata$FL), max(clim.bfdata$FL), length.out = 100) #create sequence within FL range
pred.y.mod1 = predict(clim.bf.mod1, data.frame(FL = seq.x)) #based on FL sequence, estimate ABF for model 1
pred.y.mod2 = predict(clim.bf.mod2, data.frame(FL = seq.x, FL2 = seq.x^2)) #estimate for model 2
pred.y.mod3 = predict(clim.bf.mod3, data.frame(FL = seq.x, FL2 = seq.x^2, FL3 = seq.x^3)) #estimate for model 3
pred.y.mod4 = predict(clim.bf.mod4, data.frame(FL = seq.x, FL2 = seq.x^2, FL3 = seq.x^3, FL4 = seq.x^4)) #estimate for model 4

plot(ABF ~ FL, data = clim.bfdata, pch = 16)
lines(seq.x, pred.y.mod1, lty=2, lwd=2, col="red")
lines(seq.x, pred.y.mod2, lty=2, lwd=2, col="orange")
lines(seq.x, pred.y.mod3, lty=1, lwd=3, col="blue")
lines(seq.x, pred.y.mod4, lty=2, lwd=2, col="green")



#Extract coeffs from clim.bf.model3 and create fxn to find where 2nd derivative = 0

clim.bf.int<- as.numeric(clim.bf.mod3$coefficients[1]) #intercept
clim.bf.FL<- as.numeric(clim.bf.mod3$coefficients[2]) #FL
clim.bf.FL2<- as.numeric(clim.bf.mod3$coefficients[3]) #FL^2
clim.bf.FL3<- as.numeric(clim.bf.mod3$coefficients[4]) #FL^3

clim.bf.fxn<- function(x) clim.bf.FL3*x^3 + clim.bf.FL2*x^2 + clim.bf.FL*x + clim.bf.int

#plot
ggplot(data = clim.bfdata, aes(x = FL)) + geom_point(aes(y = ABF), size = 3) + stat_function(fun = clim.bf.fxn, size = 1, color = "blue") + theme_bw()





## S. tiburo ##

stib.bfdata<- read.csv("BH BF Values.csv", header=TRUE, sep=",")


#create models
stib.bf.mod1<- lm(ABF ~ FL, data = stib.bfdata)
stib.bf.mod2<- lm(ABF ~ poly(FL, 2, raw = T), data = stib.bfdata) #quadratic term
stib.bf.mod3<- lm(ABF ~ poly(FL, 3, raw = T), data = stib.bfdata) #cubic term
stib.bf.mod4<- lm(ABF ~ poly(FL, 4, raw = T), data = stib.bfdata) #quartic term

stib.bf.mods<- list("mod1" = stib.bf.mod1, "mod2" = stib.bf.mod2, "mod3" = stib.bf.mod3, "mod4" = stib.bf.mod4)

#test for differences among models by SS (anova fxn) and AICc
anova(stib.bf.mod1, stib.bf.mod2, stib.bf.mod3, stib.bf.mod4) #model 2 > model 1; model 2 > model 3; model 4 > model 3
anova(stib.bf.mod2, stib.bf.mod4) #model 4 is best
aictab(cand.set = stib.bf.mods, second.ord = T) #model 4 is best

#Investigate the final model (Model 4)
summary(stib.bf.mod4) #provides coefficients for model, p-value for overall model, and R2


#Plot polynomial regression from each model for comparison

seq.x = seq(min(stib.bfdata$FL), max(stib.bfdata$FL), length.out = 100) #create sequence within FL range
pred.y.mod1 = predict(stib.bf.mod1, data.frame(FL = seq.x)) #based on FL sequence, estimate ABF for model 1
pred.y.mod2 = predict(stib.bf.mod2, data.frame(FL = seq.x, FL2 = seq.x^2)) #estimate for model 2
pred.y.mod3 = predict(stib.bf.mod3, data.frame(FL = seq.x, FL2 = seq.x^2, FL3 = seq.x^3)) #estimate for model 3
pred.y.mod4 = predict(stib.bf.mod4, data.frame(FL = seq.x, FL2 = seq.x^2, FL3 = seq.x^3, FL4 = seq.x^4)) #estimate for model 4

plot(ABF ~ FL, data = stib.bfdata, pch = 16)
lines(seq.x, pred.y.mod1, lty=2, lwd=2, col="red")
lines(seq.x, pred.y.mod2, lty=2, lwd=2, col="orange")
lines(seq.x, pred.y.mod3, lty=2, lwd=2, col="green")
lines(seq.x, pred.y.mod4, lty=1, lwd=3, col="blue")



#Extract coeffs from stib.bf.model3 and create fxn to find where 2nd derivative = 0

stib.bf.int<- as.numeric(stib.bf.mod4$coefficients[1]) #intercept
stib.bf.FL<- as.numeric(stib.bf.mod4$coefficients[2]) #FL
stib.bf.FL2<- as.numeric(stib.bf.mod4$coefficients[3]) #FL^2
stib.bf.FL3<- as.numeric(stib.bf.mod4$coefficients[4]) #FL^3
stib.bf.FL4<- as.numeric(stib.bf.mod4$coefficients[5]) #FL^4

stib.bf.fxn<- function(x) stib.bf.FL4*x^4 + stib.bf.FL3*x^3 + stib.bf.FL2*x^2 + stib.bf.FL*x + stib.bf.int

#plot
ggplot(data = stib.bfdata, aes(x = FL)) + geom_point(aes(y = ABF), size = 3) + stat_function(fun = stib.bf.fxn, size = 1, color = "blue") + theme_bw()










#### Find Where 2nd Deriv = 0 ####

## C. leucas ##

cleu.2d<- Deriv(cleu.bf.fxn, "x", nderiv = 2) #calculates 2nd derivative of fxn

#find root
uniroot(function(x) cleu.2d(x), c(min(cleu.bfdata$FL),max(cleu.bfdata$FL))) #120.8816 cm FL

#check by plugging back in to fxn
cleu.2d(120.8816) #checks out

#show on plot
cleu.InfPt.plot<- ggplot(data = cleu.bfdata, aes(x = FL)) + geom_point(aes(y = ABF), size = 3, color = "steelblue4") + stat_function(fun = cleu.bf.fxn, size = 1, color = "black") + theme_bw() + geom_vline(xintercept = 120.8816, color = "red", lwd = 1) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), plot.margin = unit(c(0.5,0.5,0,0.25), "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) + labs(x = "", y = "")





## C. limbatus ##

clim.2d<- Deriv(clim.bf.fxn, "x", nderiv = 2) #calculates 2nd derivative of fxn

#find root
uniroot(function(x) clim.2d(x), c(min(clim.bfdata$FL),max(clim.bfdata$FL))) #99.23374 cm FL

#check by plugging back in to fxn
clim.2d(99.23374) #checks out

#show on plot
clim.InfPt.plot<- ggplot(data = clim.bfdata, aes(x = FL)) + geom_point(aes(y = ABF), size = 3, color = "firebrick") + stat_function(fun = clim.bf.fxn, size = 1, color = "black") + theme_bw() + geom_vline(xintercept = 99.23374, color = "red", lwd=1) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), plot.margin = unit(c(0.5,0.5,0,0.25), "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) + labs(x = "", y = "ABF (N)")





## S. tiburo ##

stib.2d<- Deriv(stib.bf.fxn, "x", nderiv = 2) #calculates 2nd derivative of fxn

#find root (search only between 65 and 90 cm FL)
uniroot(function(x) stib.2d(x), c(65,90)) #82.14386 cm FL

#check by plugging back in to fxn
stib.2d(82.14386) #checks out

#show on plot
stib.InfPt.plot<- ggplot(data = stib.bfdata, aes(x = FL)) + geom_point(aes(y = ABF), size = 3, color = "darkseagreen4") + stat_function(fun = stib.bf.fxn, size = 1, color = "black") + theme_bw() + geom_vline(xintercept = 82.14386, color = "red", lwd=1) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), panel.grid = element_blank(), plot.margin = unit(c(0.5,0.5,0,0.25), "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) + labs(x = "FL (cm)", y = "")


#Plot all together

bf_InfPt.plot<- plot_grid(cleu.InfPt.plot, clim.InfPt.plot, stib.InfPt.plot, labels = "", nrow = 3, align = "v")

save_plot("ABF Inflection Points.tiff", bf_InfPt.plot, nrow = 3, base_aspect_ratio = 2)











#### If all else fails, just run RMA (SMA) regression on log-transformed data to measure scaling relationships ####

library(smatr)

## C. leucas ##


#Small C. leucas

cleu.small.sma<- sma(ABF ~ FL, data = cleu.small, log = "xy", slope.test = 2)
plot(cleu.small.sma)
cleu.small.sma #good fit (R2 = 0.81) and slope is > 2 (95% CI: 2.43,3.70)

#check resids and normality
plot(cleu.small.sma, which = "residual") #looks good
plot(cleu.small.sma, which = "qq") #looks good



#Large C. leucas

cleu.large.sma<- sma(ABF ~ FL, data = cleu.large, log = "xy", slope.test = 2)
plot(cleu.large.sma)
cleu.large.sma #okay fit (R2 = 0.89) and slope is not different from 2 (95% CI: 0.65,3.76)

#check resids and normality
plot(cleu.large.sma, which = "residual") #can't say w only 4 points
plot(cleu.large.sma, which = "qq") #single outlier; see how robust model looks (weights values by using Huber's M)


cleu.large.sma2<- sma(ABF ~ FL, data = cleu.large, log = "xy", slope.test = 2, robust = T)
plot(cleu.large.sma2) #doesn't look any different
cleu.large.sma2 #essentially no different from original model



### Plot

#create fxns for iso and allom
cleu.small.iso.fun<- function(x) (10^-1.621801)*(x^2)
cleu.small.allom.fun<- function(x) (10^-3.524023)*(x^2.997809)

cleu.large.iso.fun<- function(x) (10^-1.470752)*(x^2)
cleu.large.allom.fun<- function(x) (10^-0.5150573)*(x^1.5648202)


ggplot(data = cleu.bfdata, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = cleu.small.iso.fun, size = 1, linetype = 1, color = "blue", xlim = c(min(cleu.small$FL), max(cleu.small$FL))) + stat_function(fun = cleu.small.allom.fun, size = 1, linetype = 2, color = "blue", xlim = c(min(cleu.small$FL), max(cleu.small$FL))) + stat_function(fun = cleu.large.iso.fun, size = 1, linetype = 1, color = "red", xlim = c(min(cleu.large$FL), max(cleu.large$FL))) + stat_function(fun = cleu.large.allom.fun, size = 1, linetype = 2, color = "red", xlim = c(min(cleu.large$FL), max(cleu.large$FL))) + theme_bw()









## C. limbatus ##

clim.small<- clim.bfdata %>% filter(FL < 99.23)
clim.large<- clim.bfdata %>% filter(FL > 99.23)


#Small C. limbatus

clim.small.sma<- sma(ABF ~ FL, data = clim.small, log = "xy", slope.test = 2)
plot(clim.small.sma)
clim.small.sma #decent fit (R2 = 0.78) and slope is not different from 2 (95% CI: 1.78,3.19)

#check resids and normality
plot(clim.small.sma, which = "residual") #looks okay
plot(clim.small.sma, which = "qq") #looks okay



#Large C. limbatus

clim.large.sma<- sma(ABF ~ FL, data = clim.large, log = "xy", slope.test = 2)
plot(clim.large.sma)
clim.large.sma #alright fit (R2 = 0.61) and slope is not different from 2 (95% CI: 1.69,3.39)

#check resids and normality
plot(clim.large.sma, which = "residual") #looks good
plot(clim.large.sma, which = "qq") #single outlier; see how robust model looks (weights values by using Huber's M)


clim.large.sma2<- sma(ABF ~ FL, data = clim.large, log = "xy", slope.test = 2, robust = T)
plot(clim.large.sma2) #minute difference in elevation
clim.large.sma2 #essentially no different from original model



### Plot

#create fxns for iso and allom
clim.small.iso.fun<- function(x) (10^-1.928375)*(x^2)
clim.small.allom.fun<- function(x) (10^-2.656715)*(x^2.384287)

clim.large.iso.fun<- function(x) (10^-1.74957)*(x^2)
clim.large.allom.fun<- function(x) (10^-2.559554)*(x^2.391851)


ggplot(data = clim.bfdata, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = clim.small.iso.fun, size = 1, linetype = 1, color = "blue", xlim = c(min(clim.small$FL), max(clim.small$FL))) + stat_function(fun = clim.small.allom.fun, size = 1, linetype = 2, color = "blue", xlim = c(min(clim.small$FL), max(clim.small$FL))) + stat_function(fun = clim.large.iso.fun, size = 1, linetype = 1, color = "red", xlim = c(min(clim.large$FL), max(clim.large$FL))) + stat_function(fun = clim.large.allom.fun, size = 1, linetype = 2, color = "red", xlim = c(min(clim.large$FL), max(clim.large$FL))) + theme_bw()










## S. tiburo ##

stib.small<- stib.bfdata %>% filter(FL < 82.14)
stib.large<- stib.bfdata %>% filter(FL > 82.14)


#Small S. tiburo

stib.small.sma<- sma(ABF ~ FL, data = stib.small, log = "xy", slope.test = 2)
plot(stib.small.sma)
stib.small.sma #not great fit (R2 = 0.48) and slope > 2 (95% CI: 2.32,4.41)

#check resids and normality
plot(stib.small.sma, which = "residual") #looks okay
plot(stib.small.sma, which = "qq") #looks okay; try robust model to see if improvement



stib.small.sma2<- sma(ABF ~ FL, data = stib.small, log = "xy", slope.test = 2, robust = T)
plot(stib.small.sma2) #minute increase in slope
stib.small.sma2 #essentially no different from original model



#Large S. tiburo

stib.large.sma<- sma(ABF ~ FL, data = stib.large, log = "xy", slope.test = 2)
plot(stib.large.sma)
stib.large.sma #poor fit (R2 = 0.05) and slope is not different from 2 (95% CI: 0.49,7.27)

#check resids and normality
plot(stib.large.sma, which = "residual") #looks alright
plot(stib.large.sma, which = "qq") #looks okay



### Plot

#create fxns for iso and allom
stib.small.iso.fun<- function(x) (10^-2.314053)*(x^2)
stib.small.allom.fun<- function(x) (10^-4.484093)*(x^3.198358)

stib.large.iso.fun<- function(x) (10^-2.045178)*(x^2)
stib.large.allom.fun<- function(x) (10^-1.830639)*(x^1.890224)


ggplot(data = stib.bfdata, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = stib.small.iso.fun, size = 1, linetype = 1, color = "blue", xlim = c(min(stib.small$FL), max(stib.small$FL))) + stat_function(fun = stib.small.allom.fun, size = 1, linetype = 2, color = "blue", xlim = c(min(stib.small$FL), max(stib.small$FL))) + stat_function(fun = stib.large.iso.fun, size = 1, linetype = 1, color = "red", xlim = c(min(stib.large$FL), max(stib.large$FL))) + stat_function(fun = stib.large.allom.fun, size = 1, linetype = 2, color = "red", xlim = c(min(stib.large$FL), max(stib.large$FL))) + theme_bw()












#### Try expanding range of points analyzed by expanding beyond inflection +/- 5 cm FL for C. leucas and C. limbatus and +/- 1 cm FL for S. tiburo ####

## Start by just analyzing with 'smatr' for now


library(smatr)

## C. leucas ##

cleu.small2<- cleu.bfdata %>% filter(FL < (120.88 + 5))
cleu.large2<- cleu.bfdata %>% filter(FL > (120.88 - 5))

ggplot(data = cleu.bfdata, aes(FL)) + geom_point(data = cleu.small2, aes(y = ABF), size = 3, color = "blue", alpha = 0.8) + geom_point(data = cleu.large2, aes(y = ABF), size = 3, color = "red", alpha = 0.6) + theme_bw() + geom_vline(xintercept = 120.88)



#Small C. leucas

cleu.small.sma2<- sma(ABF ~ FL, data = cleu.small2, log = "xy", slope.test = 2)
plot(cleu.small.sma2)
cleu.small.sma2 #good fit (R2 = 0.81) and slope is > 2 (95% CI: 2.43,3.70)

#check resids and normality
plot(cleu.small.sma2, which = "residual") #looks good
plot(cleu.small.sma2, which = "qq") #looks good



#Large C. leucas

cleu.large.sma2<- sma(ABF ~ FL, data = cleu.large2, log = "xy", slope.test = 2)
plot(cleu.large.sma2)
cleu.large.sma2 #good fit (R2 = 0.95) and slope is not different from 2 (95% CI: 1.16,2.50)

#check resids and normality
plot(cleu.large.sma2, which = "residual") #looks good
plot(cleu.large.sma2, which = "qq") #single outlier; try w robust model (Huber's M)



cleu.large.sma3<- sma(ABF ~ FL, data = cleu.large2, log = "xy", slope.test = 2, robust = T)
plot(cleu.large.sma3)
cleu.large.sma3 #good fit (R2 = 0.95) and slope is not different from 2 (95% CI: 1.33,2.12)

#check resids and normality
plot(cleu.large.sma3, which = "residual") #looks good
plot(cleu.large.sma3, which = "qq") #still has outlier, but likely better estimate



### Plot

#create fxns for iso and allom
cleu.small.iso.fun2<- function(x) (10^-1.621801)*(x^2)
cleu.small.allom.fun2<- function(x) (10^-3.524023)*(x^2.997809)

cleu.large.iso.fun2<- function(x) (10^-1.465986)*(x^2)
cleu.large.allom.fun2<- function(x) (10^-0.77765161)*(x^1.680699)


ggplot(data = cleu.bfdata, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = cleu.small.iso.fun2, size = 1, linetype = 1, color = "blue", xlim = c(min(cleu.small2$FL), max(cleu.small2$FL))) + stat_function(fun = cleu.small.allom.fun2, size = 1, linetype = 2, color = "blue", xlim = c(min(cleu.small2$FL), max(cleu.small2$FL))) + stat_function(fun = cleu.large.iso.fun2, size = 1, linetype = 1, color = "red", xlim = c(min(cleu.large2$FL), max(cleu.large2$FL))) + stat_function(fun = cleu.large.allom.fun2, size = 1, linetype = 2, color = "red", xlim = c(min(cleu.large2$FL), max(cleu.large2$FL))) + theme_bw()









## C. limbatus ##

clim.small2<- clim.bfdata %>% filter(FL < 99.23 + 5)
clim.large2<- clim.bfdata %>% filter(FL > 99.23 - 5)


ggplot(data = clim.bfdata, aes(FL)) + geom_point(data = clim.small2, aes(y = ABF), size = 3, color = "blue", alpha = 0.8) + geom_point(data = clim.large2, aes(y = ABF), size = 3, color = "red", alpha = 0.6) + theme_bw() + geom_vline(xintercept = 99.23)


#Small C. limbatus

clim.small.sma2<- sma(ABF ~ FL, data = clim.small2, log = "xy", slope.test = 2)
plot(clim.small.sma2)
clim.small.sma2 #decent fit (R2 = 0.80) and slope > 2 (95% CI: 2.06,3.43)

#check resids and normality
plot(clim.small.sma2, which = "residual") #looks okay
plot(clim.small.sma2, which = "qq") #a little off at the lower tail; try robust model



clim.small.sma3<- sma(ABF ~ FL, data = clim.small2, log = "xy", slope.test = 2, robust = T)
plot(clim.small.sma3) #slight increase in slope
clim.small.sma3 #decent fit (R2 = 0.80) and slope > 2 (95% CI: 2.15, 3.90)

#check resids and normality
plot(clim.small.sma3, which = "residual") #looks okay
plot(clim.small.sma3, which = "qq") #slight improvement; use this model




#Large C. limbatus

clim.large.sma2<- sma(ABF ~ FL, data = clim.large2, log = "xy", slope.test = 2)
plot(clim.large.sma2)
clim.large.sma2 #decent fit (R2 = 0.73) and slope > 2 (95% CI: 2.01,3.45)

#check resids and normality
plot(clim.large.sma2, which = "residual") #looks good
plot(clim.large.sma2, which = "qq") #single outlier on lower tail; try robust method



clim.large.sma3<- sma(ABF ~ FL, data = clim.large2, log = "xy", slope.test = 2, robust = T)
plot(clim.large.sma3)
clim.large.sma3 #decent fit (R2 = 0.73) and slope > 2 (95% CI: 2.004,3.27)

#check resids and normality
plot(clim.large.sma3, which = "residual") #looks good
plot(clim.large.sma3, which = "qq") #no improvement; stick with original model


### Plot

#create fxns for iso and allom
clim.small.iso.fun2<- function(x) (10^-1.904921)*(x^2)
clim.small.allom.fun2<- function(x) (10^-3.617967)*(x^2.895455)

clim.large.iso.fun2<- function(x) (10^-1.759632)*(x^2)
clim.large.allom.fun2<- function(x) (10^-3.066041)*(x^2.635060)


ggplot(data = clim.bfdata, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = clim.small.iso.fun2, size = 1, linetype = 1, color = "blue", xlim = c(min(clim.small2$FL), max(clim.small2$FL))) + stat_function(fun = clim.small.allom.fun2, size = 1, linetype = 2, color = "blue", xlim = c(min(clim.small2$FL), max(clim.small2$FL))) + stat_function(fun = clim.large.iso.fun2, size = 1, linetype = 1, color = "red", xlim = c(min(clim.large2$FL), max(clim.large2$FL))) + stat_function(fun = clim.large.allom.fun2, size = 1, linetype = 2, color = "red", xlim = c(min(clim.large2$FL), max(clim.large2$FL))) + theme_bw()










## S. tiburo ##

stib.small2<- stib.bfdata %>% filter(FL < 82.14 + 1)
stib.large2<- stib.bfdata %>% filter(FL > 82.14 - 1)

ggplot(data = stib.bfdata, aes(FL)) + geom_point(data = stib.small2, aes(y = ABF), size = 3, color = "blue", alpha = 0.8) + geom_point(data = stib.large2, aes(y = ABF), size = 3, color = "red", alpha = 0.6) + theme_bw() + geom_vline(xintercept = 82.14)


#Small S. tiburo

stib.small.sma2<- sma(ABF ~ FL, data = stib.small2, log = "xy", slope.test = 2)
plot(stib.small.sma2)
stib.small.sma2 #not great fit (R2 = 0.48) and slope > 2 (95% CI: 2.32,4.41)

#check resids and normality
plot(stib.small.sma2, which = "residual") #looks okay
plot(stib.small.sma2, which = "qq") #looks good




#Large S. tiburo

stib.large.sma2<- sma(ABF ~ FL, data = stib.large2, log = "xy", slope.test = 2)
plot(stib.large.sma2)
stib.large.sma2 #poor fit (R2 = 0.08) and slope not different than 2 (95% CI: 0.54,4.83)

#check resids and normality
plot(stib.large.sma2, which = "residual") #looks alright
plot(stib.large.sma2, which = "qq") #looks okay



### Plot

#create fxns for iso and allom
stib.small.iso.fun2<- function(x) (10^-2.314053)*(x^2)
stib.small.allom.fun2<- function(x) (10^-4.484093)*(x^3.198358)

stib.large.iso.fun2<- function(x) (10^-2.034453)*(x^2)
stib.large.allom.fun2<- function(x) (10^-1.277219)*(x^1.6110855)


ggplot(data = stib.bfdata, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = stib.small.iso.fun2, size = 1, linetype = 1, color = "blue", xlim = c(min(stib.small2$FL), max(stib.small2$FL))) + stat_function(fun = stib.small.allom.fun2, size = 1, linetype = 2, color = "blue", xlim = c(min(stib.small2$FL), max(stib.small2$FL))) + stat_function(fun = stib.large.iso.fun2, size = 1, linetype = 1, color = "red", xlim = c(min(stib.large2$FL), max(stib.large2$FL))) + stat_function(fun = stib.large.allom.fun2, size = 1, linetype = 2, color = "red", xlim = c(min(stib.large2$FL), max(stib.large2$FL))) + theme_bw()













#### Calculate Scaling Relationships for Each Species via Bayesian Regression ####

## C. leucas ##

set.seed(123)

#separate DFs based upon inflection point (120.88 cm)

## Per BSBs notes, a high corr between two explanatory vars (e.g. a and b) could drastically increase the chain length needed and result in autocorr; try centering the FL


##Data
ABF<-cleu.small2$ABF
FL<- cleu.small2$FL
log.FL<- log(FL)
log.FL.c<- log.FL - mean(log.FL)
n<- nrow(cleu.small2)
cleu.small2$log.FL.c<- (log(cleu.small2$FL) - mean(log(cleu.small2$FL))) #add in DF for R2 calc later

cleu.small.data<- list("ABF","log.FL.c","n")

##Model

cleu.small.model.c<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log(a) + b*log.FL.c[k]
    bf[k]<- exp(log.bf[k])
    ABF[k] ~ dnorm(bf[k], tau)
  }
  #Priors
  a ~ dunif(0,50)
  b ~ dunif(0.01,15)
  sigma ~ dunif(0,112*2)  #Sigma chosen by finding (max(ABF)-min(ABF))/4 as center of positive distribution
  tau<- 1/sigma/sigma
}

##Params

cleu.small.params<- c("a","b")

##Inits

cleu.small.inits<- function() {
  list("a" = runif(1,0,50), "b" = runif(1,0.01,15), "sigma" = runif(1,0,112*2))
}

##Run model

cleu.small.fit.c<-jags(data=cleu.small.data, 
                       inits=cleu.small.inits, 
                       parameters.to.save=cleu.small.params,
                       model.file=cleu.small.model.c,
                       n.iter=10000, n.burnin=1, n.thin=1, n.chains=1)
cleu.small.mcmc.c<- as.mcmc(cleu.small.fit.c)

traceplot(cleu.small.fit.c) #check traceplots for a, b, and deviance
autocorr.plot(cleu.small.mcmc.c) #autocorr high
crosscorr(cleu.small.mcmc.c); crosscorr.plot(cleu.small.mcmc.c) #high neg corr between 'a' and 'b'
raftery.diag(cleu.small.mcmc.c) #check lower quantile for 95% probability
raftery.diag(cleu.small.mcmc.c, q=0.975) #check upper quantile for 95% probability

##traceplot for 'a' show it's hitting ceiling, so increase max val in distrib to 500. For deviance, raftery suggests using a burn-in of 98 and dependence factor of 27 (for thinning); re-run with new values








##Model 2

cleu.small.model.c2<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log(a) + b*log.FL.c[k]
    bf[k]<- exp(log.bf[k])
    ABF[k] ~ dnorm(bf[k], tau)
  }
  #Priors
  a ~ dunif(0,500) #changed from 50 to 500 to raise ceiling
  b ~ dunif(0.01,15)
  sigma ~ dunif(0,112*2)
  tau<- 1/sigma/sigma
  # Posterior predictive simulations 
  for (i in 1:n) {
    abfPred[i]~dnorm(bf[i], tau)
  }
}


##Params

cleu.small.params2<- c("a","b","bf","abfPred")


##Inits

cleu.small.inits2<- function() {
  #for 'a', 50 changed to 500 to raise ceiling
  list("a" = runif(1,0,500), "b" = runif(1,0.01,15), "sigma" = runif(1,0,112*2))
}

##Run model

cleu.small.fit.c2<-jags(data=cleu.small.data, 
                        inits=cleu.small.inits2, 
                        parameters.to.save=cleu.small.params2,
                        model.file=cleu.small.model.c2,
                        n.iter=10000*27+98, n.burnin=98, n.thin=27, n.chains=1)
cleu.small.mcmc.c2<- as.mcmc(cleu.small.fit.c2)

traceplot(cleu.small.fit.c2) #mixing looks good
autocorr.plot(cleu.small.mcmc.c2) #looks good
crosscorr(cleu.small.mcmc.c2); crosscorr.plot(cleu.small.mcmc.c2) #high neg corr between 'a' and 'b'



##Now run model with 3 chains and evaulate diagnostics for convergence



cleu.small.fit.c3<-jags(data=cleu.small.data, 
                        inits=cleu.small.inits2, 
                        parameters.to.save=cleu.small.params2,
                        model.file=cleu.small.model.c2,
                        n.iter=10000*27+98, n.burnin=98, n.thin=27, n.chains=3)
cleu.small.mcmc.c3<- as.mcmc(cleu.small.fit.c3)


plot(cleu.small.mcmc.c3[,c("a","b")]) #mixing looks good
acfplot(cleu.small.mcmc.c3[,c("a","b")]) #looks good
crosscorr(cleu.small.mcmc.c3[,c("a","b")]); crosscorr.plot(cleu.small.mcmc.c3[,c("a","b")]) #high neg corr between 'a' and 'b'

#Run diagnostic tests
gelman.diag(cleu.small.mcmc.c3[,c("a","b")]); gelman.plot(cleu.small.mcmc.c3[,c("a","b")]) #looks very good
geweke.diag(cleu.small.mcmc.c3[,c("a","b")]); geweke.plot(cleu.small.mcmc.c3[,c("a","b")]) #appears to show convergence



## Posterior Predictive Check (of residuals and distributions)

simulations = cleu.small.fit.c3$BUGSoutput$sims.list$abfPred
pred = apply(cleu.small.fit.c3$BUGSoutput$sims.list$bf, 2, median)
dim(simulations)
sim = createDHARMa(simulatedResponse = t(simulations), observedResponse = cleu.small2$ABF, fittedPredictedResponse = pred, integerResponse = F)
plot(sim) #normality and HOV check out
ggplot() + geom_density(data = NULL, aes(x = as.vector(simulations), fill = "Model"), alpha = 0.5) + geom_density(data = cleu.small2, aes(x = ABF, fill = "Obs"), alpha = 0.5) #matches up well



##Get summary data of power function and plot

cleu.small.fit.c3
summary(cleu.small.mcmc.c3)

#create new power fxn

a<- as.numeric(cleu.small.fit.c3$BUGSoutput$mean$a)
b<- as.numeric(cleu.small.fit.c3$BUGSoutput$mean$b)

cleu.small.fxn<- function(x) a*x^b

ggplot(data = cleu.small2, aes(x = exp(log.FL.c), y = ABF)) + geom_point(size = 3) + stat_function(fun = cleu.small.fxn, size = 1, color = "blue") + theme_bw()

modifier.for.a <- exp(-b*mean(log(cleu.small2$FL)))
a*modifier.for.a # 0.0001431931
cleu.small.fxn.alt <- function(x) a*modifier.for.a*x^b

ggplot(data = cleu.small2, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = cleu.small.fxn.alt, size = 1, color = "blue") + theme_bw()


#Bayesian R2 (via Gelman et al 2017)

mcmc.cleu.small.fit <- cleu.small.fit.c3$BUGSoutput$sims.matrix
Xmat.cleu.small.fit = model.matrix(~log.FL.c, cleu.small2)
coefs.cleu.small.fit = mcmc.cleu.small.fit[, c("a", "b")]
fit_cleu.small.fit = coefs.cleu.small.fit %*% t(Xmat.cleu.small.fit)
resid_cleu.small.fit = sweep(fit_cleu.small.fit, 2, log(cleu.small2$ABF), "-")
var_f_cleu.small = apply(fit_cleu.small.fit, 1, var)
var_e_cleu.small = apply(resid_cleu.small.fit, 1, var)
R2_cleu.small = var_f_cleu.small/(var_f_cleu.small + var_e_cleu.small)
tidyMCMC(as.mcmc(R2_cleu.small), conf.int = TRUE, conf.method = "HPDinterval")





### Cleu large model


##Data
ABF<-cleu.large2$ABF
FL<- cleu.large2$FL
log.FL<- log(FL)
log.FL.c<- log.FL - mean(log.FL)
n<- nrow(cleu.large2)
cleu.large2$log.FL.c<- (log(cleu.large2$FL) - mean(log(cleu.large2$FL)))

cleu.large.data<- list("ABF","log.FL.c","n")

##Model (centered FL)

cleu.large.model.c<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log(a) + b*log.FL.c[k]
    bf[k]<- exp(log.bf[k])
    ABF[k] ~ dnorm(bf[k], tau)
  }
  #Priors
  a ~ dunif(0,500)
  b ~ dunif(0.01,15)
  sigma ~ dunif(0,113*2)
  tau<- 1/sigma/sigma
}


##Params

cleu.large.params<- c("a","b")


##Inits

cleu.large.inits<- function() {
  #for 'a', 50 changed to 500 to raise ceiling
  list("a" = runif(1,0,500), "b" = runif(1,0.01,15), "sigma" = runif(1,0,113*2))
}

##Run model

cleu.large.fit.c<-jags(data=cleu.large.data, 
                       inits=cleu.large.inits, 
                       parameters.to.save=cleu.large.params,
                       model.file=cleu.large.model.c,
                       n.iter=10000, n.burnin=1, n.thin=1, n.chains=1)
cleu.large.mcmc.c<- as.mcmc(cleu.large.fit.c)

traceplot(cleu.large.fit.c) #'a' hit ceiling; raise to 1000
autocorr.plot(cleu.large.mcmc.c) #high autocorr
crosscorr(cleu.large.mcmc.c); crosscorr.plot(cleu.large.mcmc.c) #high neg corr between 'a' and deviance
raftery.diag(cleu.large.mcmc.c) #check lower quantile for 95% probability
raftery.diag(cleu.large.mcmc.c, q=0.975) #check upper quantile for 95% probability


#based on Raftery-Lewis, use burn-in of 152 and thinning of 42




##Model 2 (centered FL)

cleu.large.model.c2<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log(a) + b*log.FL.c[k]
    bf[k]<- exp(log.bf[k])
    ABF[k] ~ dnorm(bf[k], tau)
  }
  #Priors
  a ~ dunif(0,1000)
  b ~ dunif(0.01,15)
  sigma ~ dunif(0,113*2)
  tau<- 1/sigma/sigma
  # Posterior predictive simulations 
  for (i in 1:n) {
    abfPred[i]~dnorm(bf[i], tau)
  }
}


##Params
cleu.large.params2<- c("a","b","bf","abfPred")


##Inits

cleu.large.inits2<- function() {
  #for 'a', 50 changed to 500 to raise ceiling
  list("a" = runif(1,0,1000), "b" = runif(1,0.01,15), "sigma" = runif(1,0,113*2))
}


cleu.large.fit.c2<-jags(data=cleu.large.data, 
                        inits=cleu.large.inits2, 
                        parameters.to.save=cleu.large.params2,
                        model.file=cleu.large.model.c2,
                        n.iter=10000*42+152, n.burnin=152, n.thin=42, n.chains=1)
cleu.large.mcmc.c2<- as.mcmc(cleu.large.fit.c2)

traceplot(cleu.large.fit.c2) #mixing looks good
autocorr.plot(cleu.large.mcmc.c2) #looks good
crosscorr(cleu.large.mcmc.c2); crosscorr.plot(cleu.large.mcmc.c2) #still high for 'a' and 'b'

#now run with 3 chains



cleu.large.fit.c3<-jags(data=cleu.large.data, 
                        inits=cleu.large.inits2, 
                        parameters.to.save=cleu.large.params2,
                        model.file=cleu.large.model.c2,
                        n.iter=10000*42+152, n.burnin=152, n.thin=42, n.chains=3)
cleu.large.mcmc.c3<- as.mcmc(cleu.large.fit.c3)

plot(cleu.large.mcmc.c3[,c("a","b")]) #mixing looks good
acfplot(cleu.large.mcmc.c3[,c("a","b")]) #looks good
crosscorr(cleu.large.mcmc.c3[,c("a","b")]); crosscorr.plot(cleu.large.mcmc.c3[,c("a","b")]) #high for 'a' and 'b'

#Run diagnostic tests
gelman.diag(cleu.large.mcmc.c3[,c("a","b")]); gelman.plot(cleu.large.mcmc.c3[,c("a","b")]) #looks very good
geweke.diag(cleu.large.mcmc.c3[,c("a","b")]); geweke.plot(cleu.large.mcmc.c3[,c("a","b")]) #appears to show convergence


## Posterior Predictive Check (of residuals and distributions)

simulations = cleu.large.fit.c3$BUGSoutput$sims.list$abfPred
pred = apply(cleu.large.fit.c3$BUGSoutput$sims.list$bf, 2, median)
dim(simulations)
sim = createDHARMa(simulatedResponse = t(simulations), observedResponse = cleu.large2$ABF, fittedPredictedResponse = pred, integerResponse = F)
plot(sim) #normality and HOV check out
ggplot() + geom_density(data = NULL, aes(x = as.vector(simulations), fill = "Model"), alpha = 0.5) + geom_density(data = cleu.large2, aes(x = ABF, fill = "Obs"), alpha = 0.5) #distrib is relatively close



##Get summary data of power function and plot

cleu.large.fit.c3
summary(cleu.large.mcmc.c3)

#create new power fxn

a<- as.numeric(cleu.large.fit.c3$BUGSoutput$mean$a)
b<- as.numeric(cleu.large.fit.c3$BUGSoutput$mean$b)

cleu.large.fxn<- function(x) a*x^b

ggplot(data = cleu.large2, aes(x = exp(log.FL.c), y = ABF)) + geom_point(size = 3) + stat_function(fun = cleu.large.fxn, size = 1, color = "blue") + theme_bw()

modifier.for.a <- exp(-b*mean(log(cleu.large2$FL)))
a*modifier.for.a # 0.2630845
cleu.large.fxn.alt <- function(x) a*modifier.for.a*x^b

ggplot(data = cleu.large2, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = cleu.large.fxn.alt, size = 1, color = "blue") + theme_bw()


#Bayesian R2 (via Gelman et al 2017)

mcmc.cleu.large.fit <- cleu.large.fit.c3$BUGSoutput$sims.matrix
Xmat.cleu.large.fit = model.matrix(~log.FL.c, cleu.large2)
coefs.cleu.large.fit = mcmc.cleu.large.fit[, c("a", "b")]
fit_cleu.large.fit = coefs.cleu.large.fit %*% t(Xmat.cleu.large.fit)
resid_cleu.large.fit = sweep(fit_cleu.large.fit, 2, log(cleu.large2$ABF), "-")
var_f_cleu.large = apply(fit_cleu.large.fit, 1, var)
var_e_cleu.large = apply(resid_cleu.large.fit, 1, var)
R2_cleu.large = var_f_cleu.large/(var_f_cleu.large + var_e_cleu.large)
tidyMCMC(as.mcmc(R2_cleu.large), conf.int = TRUE, conf.method = "HPDinterval")




#create fxns for iso and allom
cleuB.small.allom.fun<- function(x) 0.0001431931*(x^3.160311)
cleuB.large.allom.fun<- function(x) 0.2630845*(x^1.591666)











## C. limbatus ##

set.seed(123)

#separate DFs based upon inflection point (99.23 cm)



##Data
ABF<-clim.small2$ABF
FL<- clim.small2$FL
log.FL<- log(FL)
log.FL.c<- log.FL - mean(log.FL)
n<- nrow(clim.small2)
clim.small2$log.FL.c<- (log(clim.small2$FL) - mean(log(clim.small2$FL)))

clim.small.data<- list("ABF","log.FL.c","n")

##Model

clim.small.model.c<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log(a) + b*log.FL.c[k]
    bf[k]<- exp(log.bf[k])
    ABF[k] ~ dnorm(bf[k], tau)
  }
  #Priors
  a ~ dunif(0,50)
  b ~ dunif(0.01,15)
  sigma ~ dunif(0,45.305*2)  #Sigma chosen by finding (max(ABF)-min(ABF))/4 as center of positive distribution
  tau<- 1/sigma/sigma
}

##Params

clim.small.params<- c("a","b")

##Inits

clim.small.inits<- function() {
  list("a" = runif(1,0,50), "b" = runif(1,0.01,15), "sigma" = runif(1,0,45.305*2))
}

##Run model

clim.small.fit.c<-jags(data=clim.small.data, 
                       inits=clim.small.inits, 
                       parameters.to.save=clim.small.params,
                       model.file=clim.small.model.c,
                       n.iter=10000, n.burnin=1, n.thin=1, n.chains=1)
clim.small.mcmc.c<- as.mcmc(clim.small.fit.c)

traceplot(clim.small.fit.c) #check traceplots for a, b, and deviance
autocorr.plot(clim.small.mcmc.c) #autocorr high
crosscorr(clim.small.mcmc.c); crosscorr.plot(clim.small.mcmc.c) #high neg corr between 'a' and 'b'
raftery.diag(clim.small.mcmc.c) #check lower quantile for 95% probability
raftery.diag(clim.small.mcmc.c, q=0.975) #check upper quantile for 95% probability

##traceplot for 'a' show it's hitting ceiling, so increase max val in distrib to 500. Raftery-Lewis suggests burn-in of 30 and thinning of 9; re-run with new values





##Model 2

clim.small.model.c2<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log(a) + b*log.FL.c[k]
    bf[k]<- exp(log.bf[k])
    ABF[k] ~ dnorm(bf[k], tau)
  }
  #Priors
  a ~ dunif(0,500) #changed from 50 to 500 to raise ceiling
  b ~ dunif(0.01,15)
  sigma ~ dunif(0,45.305*2)
  tau<- 1/sigma/sigma
  # Posterior predictive simulations 
  for (i in 1:n) {
    abfPred[i]~dnorm(bf[i], tau)
  }
}


##Params
clim.small.params2<- c("a","b","bf","abfPred")


##Inits

clim.small.inits2<- function() {
  #for 'a', 50 changed to 500 to raise ceiling
  list("a" = runif(1,0,500), "b" = runif(1,0.01,15), "sigma" = runif(1,0,45.305*2))
}

##Run model

clim.small.fit.c2<-jags(data=clim.small.data, 
                        inits=clim.small.inits2, 
                        parameters.to.save=clim.small.params2,
                        model.file=clim.small.model.c2,
                        n.iter=10000*9+30, n.burnin=30, n.thin=9, n.chains=1)
clim.small.mcmc.c2<- as.mcmc(clim.small.fit.c2)

traceplot(clim.small.fit.c2) #mixing looks good, but probably need to increase burn-in
autocorr.plot(clim.small.mcmc.c2) #still a little high
crosscorr(clim.small.mcmc.c2); crosscorr.plot(clim.small.mcmc.c2) #high neg corr between 'a' and 'b'
raftery.diag(clim.small.mcmc.c2) 
raftery.diag(clim.small.mcmc.c2, q=0.975) 


# Per R-L, use burn-in of 108 and thinning of 43


clim.small.fit.c3<-jags(data=clim.small.data, 
                        inits=clim.small.inits2, 
                        parameters.to.save=clim.small.params2,
                        model.file=clim.small.model.c2,
                        n.iter=10000*43+108, n.burnin=108, n.thin=43, n.chains=1)
clim.small.mcmc.c3<- as.mcmc(clim.small.fit.c3)

traceplot(clim.small.fit.c3[,c("a","b")]) #looks good
autocorr.plot(clim.small.mcmc.c3[,c("a","b")]) #looks good
crosscorr(clim.small.mcmc.c3[,c("a","b")]); crosscorr.plot(clim.small.mcmc.c3[,c("a","b")]) #high neg corr between 'a' and 'b'



##Now run model with 3 chains and evaulate diagnostics for convergence



clim.small.fit.c4<-jags(data=clim.small.data, 
                        inits=clim.small.inits2, 
                        parameters.to.save=clim.small.params2,
                        model.file=clim.small.model.c2,
                        n.iter=10000*43+108, n.burnin=108, n.thin=43, n.chains=3)
clim.small.mcmc.c4<- as.mcmc(clim.small.fit.c4)


plot(clim.small.mcmc.c4[,c("a","b")]) #mixing looks good
acfplot(clim.small.mcmc.c4[,c("a","b")]) #looks good for 'a' and 'b'
crosscorr(clim.small.mcmc.c4[,c("a","b")]); crosscorr.plot(clim.small.mcmc.c4[,c("a","b")]) #high neg corr between 'a' and 'b'

#Run diagnostic tests
gelman.diag(clim.small.mcmc.c4[,c("a","b")]); gelman.plot(clim.small.mcmc.c4[,c("a","b")]) #looks very good
geweke.diag(clim.small.mcmc.c4[,c("a","b")]); geweke.plot(clim.small.mcmc.c4[,c("a","b")]) #appears to show convergence


## Posterior Predictive Check (of residuals and distributions)

simulations = clim.small.fit.c4$BUGSoutput$sims.list$abfPred
pred = apply(clim.small.fit.c4$BUGSoutput$sims.list$bf, 2, median)
dim(simulations)
sim = createDHARMa(simulatedResponse = t(simulations), observedResponse = clim.small2$ABF, fittedPredictedResponse = pred, integerResponse = F)
plot(sim) #normality and HOV check out
ggplot() + geom_density(data = NULL, aes(x = as.vector(simulations), fill = "Model"), alpha = 0.5) + geom_density(data = clim.small2, aes(x = ABF, fill = "Obs"), alpha = 0.5) #distrib is relatively close



##Get summary data of power function and plot

clim.small.fit.c4
summary(clim.small.mcmc.c4)

#create new power fxn

a<- as.numeric(clim.small.fit.c4$BUGSoutput$mean$a)
b<- as.numeric(clim.small.fit.c4$BUGSoutput$mean$b)

clim.small.fxn<- function(x) a*x^b

ggplot(data = clim.small2, aes(x = exp(log.FL.c), y = ABF)) + geom_point(size = 3) + stat_function(fun = clim.small.fxn, size = 1, color = "blue") + theme_bw()

modifier.for.a <- exp(-b*mean(log(clim.small2$FL)))
a*modifier.for.a # 2.896163e-06
clim.small.fxn.alt <- function(x) a*modifier.for.a*x^b

ggplot(data = clim.small2, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = clim.small.fxn.alt, size = 1, color = "blue") + theme_bw()


#Bayesian R2 (via Gelman et al 2017)

mcmc.clim.small.fit <- clim.small.fit.c4$BUGSoutput$sims.matrix
Xmat.clim.small.fit = model.matrix(~log.FL.c, clim.small2)
coefs.clim.small.fit = mcmc.clim.small.fit[, c("a", "b")]
fit_clim.small.fit = coefs.clim.small.fit %*% t(Xmat.clim.small.fit)
resid_clim.small.fit = sweep(fit_clim.small.fit, 2, log(clim.small2$ABF), "-")
var_f_clim.small = apply(fit_clim.small.fit, 1, var)
var_e_clim.small = apply(resid_clim.small.fit, 1, var)
R2_clim.small = var_f_clim.small/(var_f_clim.small + var_e_clim.small)
tidyMCMC(as.mcmc(R2_clim.small), conf.int = TRUE, conf.method = "HPDinterval")




### Clim large model


##Data
ABF<-clim.large2$ABF
FL<- clim.large2$FL
log.FL<- log(FL)
log.FL.c<- log.FL - mean(log.FL)
n<- nrow(clim.large2)
clim.large2$log.FL.c<- (log(clim.large2$FL) - mean(log(clim.large2$FL)))

clim.large.data<- list("ABF","log.FL.c","n")

##Model (centered FL)

clim.large.model.c<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log(a) + b*log.FL.c[k]
    bf[k]<- exp(log.bf[k])
    ABF[k] ~ dnorm(bf[k], tau)
  }
  #Priors
  a ~ dunif(0,500)
  b ~ dunif(0.01,15)
  sigma ~ dunif(0,58.29*2)
  tau<- 1/sigma/sigma
  # Posterior predictive simulations 
  for (i in 1:n) {
    abfPred[i]~dnorm(bf[i], tau)
  }
}


##Params

clim.large.params<- c("a","b","bf","abfPred")


##Inits

clim.large.inits<- function() {
  #for 'a', 50 changed to 500 to raise ceiling
  list("a" = runif(1,0,500), "b" = runif(1,0.01,15), "sigma" = runif(1,0,58.29*2))
}

##Run model

clim.large.fit.c<-jags(data=clim.large.data, 
                       inits=clim.large.inits, 
                       parameters.to.save=clim.large.params,
                       model.file=clim.large.model.c,
                       n.iter=10000, n.burnin=1, n.thin=1, n.chains=1)
clim.large.mcmc.c<- as.mcmc(clim.large.fit.c)

traceplot(clim.large.fit.c) #increase burn-in and thinning
autocorr.plot(clim.large.mcmc.c) #high autocorr
crosscorr(clim.large.mcmc.c); crosscorr.plot(clim.large.mcmc.c) #high neg corr between 'a' and deviance
raftery.diag(clim.large.mcmc.c) #check lower quantile for 95% probability
raftery.diag(clim.large.mcmc.c, q=0.975) #check upper quantile for 95% probability


#based on Raftery-Lewis, use burn-in of 114 and thinning of 29





clim.large.fit.c2<-jags(data=clim.large.data, 
                        inits=clim.large.inits, 
                        parameters.to.save=clim.large.params,
                        model.file=clim.large.model.c,
                        n.iter=10000*29+114, n.burnin=114, n.thin=29, n.chains=1)
clim.large.mcmc.c2<- as.mcmc(clim.large.fit.c2)

traceplot(clim.large.fit.c2) #mixing looks good
autocorr.plot(clim.large.mcmc.c2) #looks good
crosscorr(clim.large.mcmc.c2); crosscorr.plot(clim.large.mcmc.c2) #looks good

#now run with 3 chains



clim.large.fit.c3<-jags(data=clim.large.data, 
                        inits=clim.large.inits, 
                        parameters.to.save=clim.large.params,
                        model.file=clim.large.model.c,
                        n.iter=10000*29+114, n.burnin=114, n.thin=29, n.chains=3)
clim.large.mcmc.c3<- as.mcmc(clim.large.fit.c3)

plot(clim.large.mcmc.c3[,c("a","b")]) #mixing looks good
acfplot(clim.large.mcmc.c3[,c("a","b")]) #looks good
crosscorr(clim.large.mcmc.c3[,c("a","b")]); crosscorr.plot(clim.large.mcmc.c3[,c("a","b")]) #looks good

#Run diagnostic tests
gelman.diag(clim.large.mcmc.c3[,c("a","b")]); gelman.plot(clim.large.mcmc.c3[,c("a","b")]) #looks very good
geweke.diag(clim.large.mcmc.c3[,c("a","b")]); geweke.plot(clim.large.mcmc.c3[,c("a","b")]) #appears to show convergence


## Posterior Predictive Check (of residuals and distributions)

simulations = clim.large.fit.c3$BUGSoutput$sims.list$abfPred
pred = apply(clim.large.fit.c3$BUGSoutput$sims.list$bf, 2, median)
dim(simulations)
sim = createDHARMa(simulatedResponse = t(simulations), observedResponse = clim.large2$ABF, fittedPredictedResponse = pred, integerResponse = F)
plot(sim) #normality and HOV check out
ggplot() + geom_density(data = NULL, aes(x = as.vector(simulations), fill = "Model"), alpha = 0.5) + geom_density(data = clim.large2, aes(x = ABF, fill = "Obs"), alpha = 0.5) #distrib is relatively close



##Get summary data of power function and plot

clim.large.fit.c3
summary(clim.large.mcmc.c3)

#create new power fxn

a<- as.numeric(clim.large.fit.c3$BUGSoutput$mean$a)
b<- as.numeric(clim.large.fit.c3$BUGSoutput$mean$b)

clim.large.fxn<- function(x) a*x^b

ggplot(data = clim.large2, aes(x = exp(log.FL.c), y = ABF)) + geom_point(size = 3) + stat_function(fun = clim.large.fxn, size = 1, color = "blue") + theme_bw()

modifier.for.a <- exp(-b*mean(log(clim.large2$FL)))
a*modifier.for.a # 0.02280475
clim.large.fxn.alt <- function(x) a*modifier.for.a*x^b

ggplot(data = clim.large2, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = clim.large.fxn.alt, size = 1, color = "blue") + theme_bw()


#Bayesian R2 (via Gelman et al 2017)

mcmc.clim.large.fit <- clim.large.fit.c3$BUGSoutput$sims.matrix
Xmat.clim.large.fit = model.matrix(~log.FL.c, clim.large2)
coefs.clim.large.fit = mcmc.clim.large.fit[, c("a", "b")]
fit_clim.large.fit = coefs.clim.large.fit %*% t(Xmat.clim.large.fit)
resid_clim.large.fit = sweep(fit_clim.large.fit, 2, log(clim.large2$ABF), "-")
var_f_clim.large = apply(fit_clim.large.fit, 1, var)
var_e_clim.large = apply(resid_clim.large.fit, 1, var)
R2_clim.large = var_f_clim.large/(var_f_clim.large + var_e_clim.large)
tidyMCMC(as.mcmc(R2_clim.large), conf.int = TRUE, conf.method = "HPDinterval")



#create fxns for iso and allom
climB.small.allom.fun<- function(x) 2.896163e-06*(x^3.874413)
climB.large.allom.fun<- function(x) 0.02280475*(x^1.946165)













## S. tiburo ##

set.seed(123)

#separate DFs based upon inflection point (99.23 cm)



##Data
ABF<-stib.small2$ABF
FL<- stib.small2$FL
log.FL<- log(FL)
log.FL.c<- log.FL - mean(log.FL)
n<- nrow(stib.small2)
stib.small2$log.FL.c<- (log(stib.small2$FL) - mean(log(stib.small2$FL)))

stib.small.data<- list("ABF","log.FL.c","n")

##Model

stib.small.model.c<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log(a) + b*log.FL.c[k]
    bf[k]<- exp(log.bf[k])
    ABF[k] ~ dnorm(bf[k], tau)
  }
  #Priors
  a ~ dunif(0,50)
  b ~ dunif(0.01,15)
  sigma ~ dunif(0,16*2)  #Sigma chosen by finding (max(ABF)-min(ABF))/4 as center of positive distribution
  tau<- 1/sigma/sigma
  # Posterior predictive simulations 
  for (i in 1:n) {
    abfPred[i]~dnorm(bf[i], tau)
  }
}

##Params

stib.small.params<- c("a","b","bf","abfPred")

##Inits

stib.small.inits<- function() {
  list("a" = runif(1,0,50), "b" = runif(1,0.01,15), "sigma" = runif(1,0,16*2))
}

##Run model

stib.small.fit.c<-jags(data=stib.small.data, 
                       inits=stib.small.inits, 
                       parameters.to.save=stib.small.params,
                       model.file=stib.small.model.c,
                       n.iter=10000, n.burnin=1, n.thin=1, n.chains=1)
stib.small.mcmc.c<- as.mcmc(stib.small.fit.c)

traceplot(stib.small.fit.c) #check traceplots for a, b, and deviance
autocorr.plot(stib.small.mcmc.c) #autocorr high
crosscorr(stib.small.mcmc.c); crosscorr.plot(stib.small.mcmc.c) #high neg corr between 'a' and 'b'
raftery.diag(stib.small.mcmc.c) #check lower quantile for 95% probability
raftery.diag(stib.small.mcmc.c, q=0.975) #check upper quantile for 95% probability

##Raftery-Lewis suggests burn-in of 156 and thinning of 42; re-run with new values





stib.small.fit.c2<-jags(data=stib.small.data, 
                        inits=stib.small.inits, 
                        parameters.to.save=stib.small.params,
                        model.file=stib.small.model.c,
                        n.iter=10000*42+156, n.burnin=156, n.thin=42, n.chains=1)
stib.small.mcmc.c2<- as.mcmc(stib.small.fit.c2)

traceplot(stib.small.fit.c2) #mixing looks good
autocorr.plot(stib.small.mcmc.c2) #looks good
crosscorr(stib.small.mcmc.c2); crosscorr.plot(stib.small.mcmc.c2) #high neg corr between 'a' and 'b'


##Now run model with 3 chains and evaulate diagnostics for convergence



stib.small.fit.c3<-jags(data=stib.small.data, 
                        inits=stib.small.inits, 
                        parameters.to.save=stib.small.params,
                        model.file=stib.small.model.c,
                        n.iter=10000*42+156, n.burnin=156, n.thin=42, n.chains=3)
stib.small.mcmc.c3<- as.mcmc(stib.small.fit.c3)


plot(stib.small.mcmc.c3[,c("a","b")]) #mixing looks good
acfplot(stib.small.mcmc.c3[,c("a","b")]) #looks good
crosscorr(stib.small.mcmc.c3[,c("a","b")]); crosscorr.plot(stib.small.mcmc.c3[,c("a","b")]) #still high between 'a' and 'b'

#Run diagnostic tests
gelman.diag(stib.small.mcmc.c3[,c("a","b")]); gelman.plot(stib.small.mcmc.c3[,c("a","b")]) #looks very good
geweke.diag(stib.small.mcmc.c3[,c("a","b")]); geweke.plot(stib.small.mcmc.c3[,c("a","b")]) #appears to show convergence


## Posterior Predictive Check (of residuals and distributions)

simulations = stib.small.fit.c3$BUGSoutput$sims.list$abfPred
pred = apply(stib.small.fit.c3$BUGSoutput$sims.list$bf, 2, median)
dim(simulations)
sim = createDHARMa(simulatedResponse = t(simulations), observedResponse = stib.small2$ABF, fittedPredictedResponse = pred, integerResponse = F)
plot(sim) #normality and HOV check out
ggplot() + geom_density(data = NULL, aes(x = as.vector(simulations), fill = "Model"), alpha = 0.5) + geom_density(data = stib.small2, aes(x = ABF, fill = "Obs"), alpha = 0.5) #distrib generally describes large peak, but this regression didn't fit as well anyway


##Get summary data of power function and plot

stib.small.fit.c3
summary(stib.small.mcmc.c3)

#create new power fxn

a<- as.numeric(stib.small.fit.c3$BUGSoutput$mean$a)
b<- as.numeric(stib.small.fit.c3$BUGSoutput$mean$b)

stib.small.fxn<- function(x) a*x^b

ggplot(data = stib.small2, aes(x = exp(log.FL.c), y = ABF)) + geom_point(size = 3) + stat_function(fun = stib.small.fxn, size = 1, color = "blue") + theme_bw()

modifier.for.a <- exp(-b*mean(log(stib.small2$FL)))
a*modifier.for.a # 2.730506e-07
stib.small.fxn.alt <- function(x) a*modifier.for.a*x^b

ggplot(data = stib.small2, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = stib.small.fxn.alt, size = 1, color = "blue") + theme_bw()


#Bayesian R2 (via Gelman et al 2017)

mcmc.stib.small.fit <- stib.small.fit.c3$BUGSoutput$sims.matrix
Xmat.stib.small.fit = model.matrix(~log.FL.c, stib.small2)
coefs.stib.small.fit = mcmc.stib.small.fit[, c("a", "b")]
fit_stib.small.fit = coefs.stib.small.fit %*% t(Xmat.stib.small.fit)
resid_stib.small.fit = sweep(fit_stib.small.fit, 2, log(stib.small2$ABF), "-")
var_f_stib.small = apply(fit_stib.small.fit, 1, var)
var_e_stib.small = apply(resid_stib.small.fit, 1, var)
R2_stib.small = var_f_stib.small/(var_f_stib.small + var_e_stib.small)
tidyMCMC(as.mcmc(R2_stib.small), conf.int = TRUE, conf.method = "HPDinterval")





### Stib large model


##Data
ABF<-stib.large2$ABF
FL<- stib.large2$FL
log.FL<- log(FL)
log.FL.c<- log.FL - mean(log.FL)
n<- nrow(stib.large2)
stib.large2$log.FL.c<- (log(stib.large2$FL) - mean(log(stib.large2$FL)))

stib.large.data<- list("ABF","log.FL.c","n")

##Model (centered FL)

stib.large.model.c<- function() {
  #Likelihood
  for (k in 1:n) {
    log.bf[k]<- log(a) + b*log.FL.c[k]
    bf[k]<- exp(log.bf[k])
    ABF[k] ~ dnorm(bf[k], tau)
  }
  #Priors
  a ~ dunif(0,500)
  b ~ dunif(0.01,15)
  sigma ~ dunif(0,6*2)
  tau<- 1/sigma/sigma
  # Posterior predictive simulations 
  for (i in 1:n) {
    abfPred[i]~dnorm(bf[i], tau)
  }
}


##Params

stib.large.params<- c("a","b","bf","abfPred")


##Inits

stib.large.inits<- function() {
  #for 'a', 50 changed to 500 to raise ceiling
  list("a" = runif(1,0,500), "b" = runif(1,0.01,15), "sigma" = runif(1,0,6*2))
}

##Run model

stib.large.fit.c<-jags(data=stib.large.data, 
                       inits=stib.large.inits, 
                       parameters.to.save=stib.large.params,
                       model.file=stib.large.model.c,
                       n.iter=10000, n.burnin=1, n.thin=1, n.chains=1)
stib.large.mcmc.c<- as.mcmc(stib.large.fit.c)

traceplot(stib.large.fit.c) #increase burn-in and thinning
autocorr.plot(stib.large.mcmc.c) #high autocorr
crosscorr(stib.large.mcmc.c); crosscorr.plot(stib.large.mcmc.c) #looks good
raftery.diag(stib.large.mcmc.c) #check lower quantile for 95% probability
raftery.diag(stib.large.mcmc.c, q=0.975) #check upper quantile for 95% probability


#based on Raftery-Lewis, use burn-in of 100 and thinning of 14





stib.large.fit.c2<-jags(data=stib.large.data, 
                        inits=stib.large.inits, 
                        parameters.to.save=stib.large.params,
                        model.file=stib.large.model.c,
                        n.iter=10000*14+100, n.burnin=100, n.thin=14, n.chains=1)
stib.large.mcmc.c2<- as.mcmc(stib.large.fit.c2)

traceplot(stib.large.fit.c2) #mixing looks good
autocorr.plot(stib.large.mcmc.c2) #looks good
crosscorr(stib.large.mcmc.c2); crosscorr.plot(stib.large.mcmc.c2) #looks good

#now run with 3 chains



stib.large.fit.c3<-jags(data=stib.large.data, 
                        inits=stib.large.inits, 
                        parameters.to.save=stib.large.params,
                        model.file=stib.large.model.c,
                        n.iter=10000*14+100, n.burnin=100, n.thin=14, n.chains=3)
stib.large.mcmc.c3<- as.mcmc(stib.large.fit.c3)

plot(stib.large.mcmc.c3[,c("a","b")]) #mixing looks good
acfplot(stib.large.mcmc.c3[,c("a","b")]) #looks good
crosscorr(stib.large.mcmc.c3[,c("a","b")]); crosscorr.plot(stib.large.mcmc.c3[,c("a","b")]) #looks good

#Run diagnostic tests
gelman.diag(stib.large.mcmc.c3[,c("a","b")]); gelman.plot(stib.large.mcmc.c3[,c("a","b")]) #looks very good
geweke.diag(stib.large.mcmc.c3[,c("a","b")]); geweke.plot(stib.large.mcmc.c3[,c("a","b")]) #appears to show convergence


## Posterior Predictive Check (of residuals and distributions)

simulations = stib.large.fit.c3$BUGSoutput$sims.list$abfPred
pred = apply(stib.large.fit.c3$BUGSoutput$sims.list$bf, 2, median)
dim(simulations)
sim = createDHARMa(simulatedResponse = t(simulations), observedResponse = stib.large2$ABF, fittedPredictedResponse = pred, integerResponse = F)
plot(sim) #normality and HOV check out
ggplot() + geom_density(data = NULL, aes(x = as.vector(simulations), fill = "Model"), alpha = 0.5) + geom_density(data = stib.large2, aes(x = ABF, fill = "Obs"), alpha = 0.5) #distrib is relatively close


##Get summary data of power function and plot

stib.large.fit.c3
summary(stib.large.mcmc.c3)

#create new power fxn

a<- as.numeric(stib.large.fit.c3$BUGSoutput$mean$a)
b<- as.numeric(stib.large.fit.c3$BUGSoutput$mean$b)

stib.large.fxn<- function(x) a*x^b

ggplot(data = stib.large2, aes(x = exp(log.FL.c), y = ABF)) + geom_point(size = 3) + stat_function(fun = stib.large.fxn, size = 1, color = "blue") + theme_bw()

modifier.for.a <- exp(-b*mean(log(stib.large2$FL)))
a*modifier.for.a # 2.903973
stib.large.fxn.alt <- function(x) a*modifier.for.a*x^b

ggplot(data = stib.large2, aes(x = FL, y = ABF)) + geom_point(size = 3) + stat_function(fun = stib.large.fxn.alt, size = 1, color = "blue") + theme_bw()


#Bayesian R2 (via Gelman et al 2017)

mcmc.stib.large.fit <- stib.large.fit.c3$BUGSoutput$sims.matrix
Xmat.stib.large.fit = model.matrix(~log.FL.c, stib.large2)
coefs.stib.large.fit = mcmc.stib.large.fit[, c("a", "b")]
fit_stib.large.fit = coefs.stib.large.fit %*% t(Xmat.stib.large.fit)
resid_stib.large.fit = sweep(fit_stib.large.fit, 2, log(stib.large2$ABF), "-")
var_f_stib.large = apply(fit_stib.large.fit, 1, var)
var_e_stib.large = apply(resid_stib.large.fit, 1, var)
R2_stib.large = var_f_stib.large/(var_f_stib.large + var_e_stib.large)
tidyMCMC(as.mcmc(R2_stib.large), conf.int = TRUE, conf.method = "HPDinterval")



#create fxns for iso and allom
stibB.small.allom.fun<- function(x) 2.730506e-07*(x^4.315665)
stibB.large.allom.fun<- function(x) 2.903973*(x^0.7176678)













#### Plot everything together ####

## Cleu
ggplot(data = cleu.bfdata, aes(x = FL, y = ABF)) + geom_point(size = 3, color = "steelblue4") + stat_function(fun = cleu.small.iso.fun2, size = 1, linetype = 1, color = "black", xlim = c(min(cleu.small2$FL), max(cleu.small2$FL))) + stat_function(fun = cleu.small.allom.fun2, size = 1, linetype = 3, color = "black", xlim = c(min(cleu.small2$FL), max(cleu.small2$FL))) + stat_function(fun = cleuB.small.allom.fun, size = 1, linetype = 2, color = "black", xlim = c(min(cleu.small2$FL), max(cleu.small2$FL))) + stat_function(fun = cleu.large.iso.fun2, size = 1, linetype = 1, color = "gray50", xlim = c(min(cleu.large2$FL), max(cleu.large2$FL))) + stat_function(fun = cleu.large.allom.fun2, size = 1, linetype = 3, color = "gray50", xlim = c(min(cleu.large2$FL), max(cleu.large2$FL))) + stat_function(fun = cleuB.large.allom.fun, size = 1, linetype = 2, color = "gray50", xlim = c(min(cleu.large2$FL), max(cleu.large2$FL))) + theme_bw() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), plot.margin = unit(c(0.5,0.5,0,0.25), "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) + labs(x = "FL (cm)", y = "ABF (N)")

ggsave("Cleu ABF scaling plot.png", width = 6, height = 4, units = "in", dpi = 600)



## Clim
ggplot(data = clim.bfdata, aes(x = FL, y = ABF)) + geom_point(size = 3, color = "firebrick") + stat_function(fun = clim.small.iso.fun2, size = 1, linetype = 1, color = "black", xlim = c(min(clim.small2$FL), max(clim.small2$FL))) + stat_function(fun = clim.small.allom.fun2, size = 1, linetype = 3, color = "black", xlim = c(min(clim.small2$FL), max(clim.small2$FL))) + stat_function(fun = climB.small.allom.fun, size = 1, linetype = 2, color = "black", xlim = c(min(clim.small2$FL), max(clim.small2$FL))) + stat_function(fun = clim.large.iso.fun2, size = 1, linetype = 1, color = "gray50", xlim = c(min(clim.large2$FL), max(clim.large2$FL))) + stat_function(fun = clim.large.allom.fun2, size = 1, linetype = 3, color = "gray50", xlim = c(min(clim.large2$FL), max(clim.large2$FL))) + stat_function(fun = climB.large.allom.fun, size = 1, linetype = 2, color = "gray50", xlim = c(min(clim.large2$FL), max(clim.large2$FL))) + theme_bw() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), plot.margin = unit(c(0.5,0.5,0,0.25), "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) + labs(x = "FL (cm)", y = "ABF (N)")

ggsave("Clim ABF scaling plot.png", width = 6, height = 4, units = "in", dpi = 600)



## Stib
ggplot(data = stib.bfdata, aes(x = FL, y = ABF)) + geom_point(size = 3, color = "darkseagreen4") + stat_function(fun = stib.small.iso.fun2, size = 1, linetype = 1, color = "black", xlim = c(min(stib.small2$FL), max(stib.small2$FL))) + stat_function(fun = stib.small.allom.fun2, size = 1, linetype = 3, color = "black", xlim = c(min(stib.small2$FL), max(stib.small2$FL))) + stat_function(fun = stibB.small.allom.fun, size = 1, linetype = 2, color = "black", xlim = c(min(stib.small2$FL), max(stib.small2$FL))) + stat_function(fun = stib.large.iso.fun2, size = 1, linetype = 1, color = "gray50", xlim = c(min(stib.large2$FL), max(stib.large2$FL))) + stat_function(fun = stib.large.allom.fun2, size = 1, linetype = 3, color = "gray50", xlim = c(min(stib.large2$FL), max(stib.large2$FL))) + stat_function(fun = stibB.large.allom.fun, size = 1, linetype = 2, color = "gray50", xlim = c(min(stib.large2$FL), max(stib.large2$FL))) + theme_bw() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), plot.margin = unit(c(0.5,0.5,0,0.25), "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) + labs(x = "FL (cm)", y = "ABF (N)")

ggsave("Stib ABF scaling plot.png", width = 6, height = 4, units = "in", dpi = 600)
