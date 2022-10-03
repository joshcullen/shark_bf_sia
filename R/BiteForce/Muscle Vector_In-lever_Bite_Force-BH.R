#For bonnetheads
#Distance between two points in 3D

setwd("~/Documents/Shark Project/")
bite <- read.table("Sharpnose Digitized Pts.csv", header=TRUE, sep=",")

#Distance from origin (x1,y1,z1) to insertion (x2,y2,z2); for muscle vectors
#Distance of Muscle (dm) = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)

#Distance from insertion (x2,y2,z2) to jaw joint (xj,yj,zj); for resolved in-lever calc
#Distance of In-Lever (dil) = sqrt((x2-xj)^2+(y2-yj)^2+(z2-zj)^2)

#Setting up coordinates
POVox<- mean(bite[1:3,2])
POVoy<- mean(bite[1:3,3])
POVoz<- mean(bite[1:3,4])
POVo<- c(POVox, POVoy, POVoz)

POVix<- mean(bite[4:6,2])
POViy<- mean(bite[4:6,3])
POViz<- mean(bite[4:6,4])
POVi<- c(POVix, POViy, POViz)

PODox<- mean(bite[7:9,2])
PODoy<- mean(bite[7:9,3])
PODoz<- mean(bite[7:9,4])
PODo<- c(PODox, PODoy, PODoz)

PODix<- mean(bite[10:12,2])
PODiy<- mean(bite[10:12,3])
PODiz<- mean(bite[10:12,4])
PODi<- c(PODix, PODiy, PODiz)

QDox<- mean(bite[13:15,2])
QDoy<- mean(bite[13:15,3])
QDoz<- mean(bite[13:15,4])
QDo<- c(QDox, QDoy, QDoz)

QDix<- mean(bite[16:18,2])
QDiy<- mean(bite[16:18,3])
QDiz<- mean(bite[16:18,4])
QDi<- c(QDix, QDiy, QDiz)

QVox<- mean(bite[19:21,2])
QVoy<- mean(bite[19:21,3])
QVoz<- mean(bite[19:21,4])
QVo<- c(QVox, QVoy, QVoz)

QVix<- mean(bite[22:24,2])
QViy<- mean(bite[22:24,3])
QViz<- mean(bite[22:24,4])
QVi<- c(QVix, QViy, QViz)

JJx<- mean(bite[25:27,2])
JJy<- mean(bite[25:27,3])
JJz<- mean(bite[25:27,4])
JJ<- c(JJx, JJy, JJz)

PBPx<- mean(bite[28:30,2])
PBPy<- mean(bite[28:30,3])
PBPz<- mean(bite[28:30,4])
PBP<- c( PBPx, PBPy, PBPz)

ABPx<- mean(bite[31:33,2])
ABPy<- mean(bite[31:33,3])
ABPz<- mean(bite[31:33,4])
ABP<- c( ABPx, ABPy, ABPz)



#Muscle Distances (dm)
POVdm<- sqrt((POVox-POVix)^2+(POVoy-POViy)^2+(POVoz-POViz)^2)
PODdm<- sqrt((PODox-PODix)^2+(PODoy-PODiy)^2+(PODoz-PODiz)^2)
QDdm<- sqrt((QDox-QDix)^2+(QDoy-QDiy)^2+(QDoz-QDiz)^2)
QVdm<- sqrt((QVox-QVix)^2+(QVoy-QViy)^2+(QVoz-QViz)^2)

#In-lever Distances (dil)
POVdil<- sqrt((JJx-POVix)^2+(JJy-POViy)^2+(JJz-POViz)^2)
PODdil<- sqrt((JJx-PODix)^2+(JJy-PODiy)^2+(JJz-PODiz)^2)
QDdil<- sqrt((JJx-QDix)^2+(JJy-QDiy)^2+(JJz-QDiz)^2)
QVdil<- sqrt((JJx-QVix)^2+(JJy-QViy)^2+(JJz-QViz)^2)



#Muscle a-CSAs
POVcsa<- 1.09
PODcsa<- 0.710
QDcsa<- 2.265
QVcsa<- 3.447
TOTALcsa<- POVcsa+PODcsa+QDcsa+QVcsa

#Individual Muscle Force Production (Po)--Bilateral 
Fpov<-POVcsa*28.9*2
Fpod<- PODcsa*28.9*2
Fqd<- QDcsa*28.9*2
Fqv<- QVcsa*28.9*2
Ftot<- Fpov+Fpod+Fqd+Fqv

#Fraction of total Po
POVf<- POVcsa/TOTALcsa
PODf<- PODcsa/TOTALcsa
QDf<- QDcsa/TOTALcsa
QVf<- QVcsa/TOTALcsa

# Resolved in-lever (weighted)
POVw<- POVf*POVdil
PODw<- PODf*PODdil
QDw<- QDf*QDdil
QVw<- QVf*QVdil

#Resolved jaw adductor in-lever
JAdil<- POVw+PODw+QDw+QVw

#Anterior Bite Point Out-lever
ABPdol<- sqrt((JJx-ABPx)^2+(JJy-ABPy)^2+(JJz-ABPz)^2)

#Posterior Bite Point Out-lever
PBPdol<- sqrt((JJx-PBPx)^2+(JJy-PBPy)^2+(JJz-PBPz)^2)

#Anterior Bite Force Production
ABF<- (Ftot*JAdil/ABPdol)

#Posterior Bite Force Production
PBF<- (Ftot*JAdil/PBPdol)




#Physiological Bite Force Estimates

#Muscle p-CSAs [=(muscle mass*cos theta)/(fiber length*muscle density)]
POVpcsa<- (2.6*1)/(4.483*1.05)
PODpcsa<- (1.2*1)/(4.483*1.05)
QDpcsa<- (3.5*1)/(4.483*1.05)  #Will likely need to alter
QVpcsa<- (6.7*1)/(4.483*1.05)
TOTALpcsa<- POVpcsa+PODpcsa+QDpcsa+QVpcsa

#Muscle Force Production (Po)--Bilateral (Physiol) 
Fpovp<-POVpcsa*28.9*2
Fpodp<- PODpcsa*28.9*2
Fqdp<- QDpcsa*28.9*2
Fqvp<- QVpcsa*28.9*2
Ftotp<- Fpovp+Fpodp+Fqdp+Fqvp

#Fraction of total Po
POVfp<- POVpcsa/TOTALpcsa
PODfp<- PODpcsa/TOTALpcsa
QDfp<- QDpcsa/TOTALpcsa
QVfp<- QVpcsa/TOTALpcsa

# Resolved in-lever (weighted)
POVwp<- POVfp*POVdil
PODwp<- PODfp*PODdil
QDwp<- QDfp*QDdil
QVwp<- QVfp*QVdil

#Resolved jaw adductor in-lever
JAdilp<- POVwp+PODwp+QDwp+QVwp


#Anterior Bite Force Production
ABFp<- (Ftotp*JAdilp/ABPdol)

#Posterior Bite Force Production
PBFp<- (Ftotp*JAdilp/PBPdol)

