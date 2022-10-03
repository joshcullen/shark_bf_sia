#### Segmented (Piece-wise) regression of Bite Force over TL ####

## need to find inflection points in model of choice (linear, polynomial, loess, etc) and then check for allometric scaling of BF

library(dplyr)
library(ggplot2)
library(inflection)


#Use BH as test example since ABF-TL relationship appears to be well-established

setwd("~/Documents/Shark Project/BF Data")
stib<- read.csv("BH BF Values.csv", header = T, sep=",")


#plot data
ggplot(stib, aes(TL,ABF)) + geom_point(size=3, color="goldenrod") + theme_bw()


### use bisection extremum distance estimator (BEDE) method to determine inflection point

inflect.stib<- bede(x = stib$TL, y = stib$ABF, index = 0)$iplast



#plot relationship with lines denoting inflection
ggplot(stib, aes(TL,ABF)) + geom_point(size=3, color="goldenrod") + theme_bw() + geom_vline(xintercept = inflect, linetype = "dashed")

#this estimated point appears to match well with the data




#### What does this look like for blacktips? ####

clim<- read.csv("BT BF Values.csv", header = T, sep=",")

#plot data
ggplot(clim, aes(TL,ABF)) + geom_point(size=3, color="indianred") + theme_bw()


### use bisection extremum distance estimator (BEDE) method to determine inflection point

inflect.clim<- bede(x = clim$TL, y = clim$ABF, index = 0)$iplast



#plot relationship with lines denoting inflection
ggplot(clim, aes(TL,ABF)) + geom_point(size=3, color="indianred") + theme_bw() + geom_vline(xintercept = inflect, linetype = "dashed")

#some high outliers for adults may be affecting the result




#### result for bull sharks? ####

cleu<- read.csv("Cleu BF Values.csv", header = T, sep=",")

#plot data
ggplot(cleu, aes(TL,ABF)) + geom_point(size=3, color="cornflowerblue") + theme_bw()


### use bisection extremum distance estimator (BEDE) method to determine inflection point

inflect.cleu<- bede(x = cleu$TL, y = cleu$ABF, index = 0)$iplast



#plot relationship with lines denoting inflection
ggplot(cleu, aes(TL,ABF)) + geom_point(size=3, color="cornflowerblue") + theme_bw() + geom_vline(xintercept = inflect, linetype = "dashed")

#inflection point looks closer to what I would guess
