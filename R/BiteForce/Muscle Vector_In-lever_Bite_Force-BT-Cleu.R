#For bull and blacktip sharks
#Distance between two points in 3D

setwd("~/Documents/Shark Project/Digitized Points/")
bite <- read.table("Cleu_Ara08_DigPts.csv", header=TRUE, sep=",")

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

QD12ox<- mean(bite[7:9,2])
QD12oy<- mean(bite[7:9,3])
QD12oz<- mean(bite[7:9,4])
QD12o<- c(QD12ox, QD12oy, QD12oz)

QD12ix<- mean(bite[10:12,2])
QD12iy<- mean(bite[10:12,3])
QD12iz<- mean(bite[10:12,4])
QD12i<- c(QD12ix, QD12iy, QD12iz)

QVox<- mean(bite[13:15,2])
QVoy<- mean(bite[13:15,3])
QVoz<- mean(bite[13:15,4])
QVo<- c(QVox, QVoy, QVoz)

QVix<- mean(bite[16:18,2])
QViy<- mean(bite[16:18,3])
QViz<- mean(bite[16:18,4])
QVi<- c(QVix, QViy, QViz)

PODox<- mean(bite[19:21,2])
PODoy<- mean(bite[19:21,3])
PODoz<- mean(bite[19:21,4])
PODo<- c(PODox, PODoy, PODoz)

PODix<- mean(bite[22:24,2])
PODiy<- mean(bite[22:24,3])
PODiz<- mean(bite[22:24,4])
PODi<- c(PODix, PODiy, PODiz)

QD3ox<- mean(bite[25:27,2])
QD3oy<- mean(bite[25:27,3])
QD3oz<- mean(bite[25:27,4])
QD3o<- c(QD3ox, QD3oy, QD3oz)

QD3ix<- mean(bite[28:30,2])
QD3iy<- mean(bite[28:30,3])
QD3iz<- mean(bite[28:30,4])
QD3i<- c(QD3ix, QD3iy, QD3iz)

QD4ox<- mean(bite[31:33,2])
QD4oy<- mean(bite[31:33,3])
QD4oz<- mean(bite[31:33,4])
QD4o<- c(QD4ox, QD4oy, QD4oz)

QD4ix<- mean(bite[34:36,2])
QD4iy<- mean(bite[34:36,3])
QD4iz<- mean(bite[34:36,4])
QD4i<- c(QD4ix, QD4iy, QD4iz)

ABPx<- mean(bite[37:39,2])
ABPy<- mean(bite[37:39,3])
ABPz<- mean(bite[37:39,4])
ABP<- c( ABPx, ABPy, ABPz)

PBPx<- mean(bite[40:42,2])
PBPy<- mean(bite[40:42,3])
PBPz<- mean(bite[40:42,4])
PBP<- c( PBPx, PBPy, PBPz)

JJx<- mean(bite[43:45,2])
JJy<- mean(bite[43:45,3])
JJz<- mean(bite[43:45,4])
JJ<- c(JJx, JJy, JJz)


#Muscle Distances (dm)
POVdm<- sqrt((POVox-POVix)^2+(POVoy-POViy)^2+(POVoz-POViz)^2)
PODdm<- sqrt((PODox-PODix)^2+(PODoy-PODiy)^2+(PODoz-PODiz)^2)
QD12dm<- sqrt((QD12ox-QD12ix)^2+(QD12oy-QD12iy)^2+(QD12oz-QD12iz)^2)
QD3dm<- sqrt((QD3ox-QD3ix)^2+(QD3oy-QD3iy)^2+(QD3oz-QD3iz)^2)
QD4dm<- sqrt((QD4ox-QD4ix)^2+(QD4oy-QD4iy)^2+(QD4oz-QD4iz)^2)
QVdm<- sqrt((QVox-QVix)^2+(QVoy-QViy)^2+(QVoz-QViz)^2)

#In-lever Distances (dil)
POVdil<- sqrt((JJx-POVix)^2+(JJy-POViy)^2+(JJz-POViz)^2)
PODdil<- sqrt((JJx-PODix)^2+(JJy-PODiy)^2+(JJz-PODiz)^2)
QD12dil<- sqrt((JJx-QD12ix)^2+(JJy-QD12iy)^2+(JJz-QD12iz)^2)
QD3dil<- sqrt((JJx-QD3ix)^2+(JJy-QD3iy)^2+(JJz-QD3iz)^2)
QD4dil<- sqrt((JJx-QD4ix)^2+(JJy-QD4iy)^2+(JJz-QD4iz)^2)
QVdil<- sqrt((JJx-QVix)^2+(JJy-QViy)^2+(JJz-QViz)^2)



#Muscle a-CSAs
POVcsa<- 2.117
PODcsa<- 0.663
QD12csa<- 0.998
QD3csa<- 1.498
QD4csa<- 0.721
QVcsa<- 5.094
TOTALcsa<- POVcsa+PODcsa+QD12csa+QD3csa+QD4csa+QVcsa

#Individual Muscle Force Production (Po)--Bilateral 
Fpov<-POVcsa*28.9*2
Fpod<- PODcsa*28.9*2
Fqd12<- QD12csa*28.9*2
Fqd3<- QD3csa*28.9*2
Fqd4<- QD4csa*14.2*2 #via citation in shark text
Fqv<- QVcsa*28.9*2
Ftot<- Fpov+Fpod+Fqd12+Fqd3+Fqd4+Fqv

#Fraction of total Po - used for in-lever calc
POVf<- Fpov/Ftot
PODf<- Fpod/Ftot
QD12f<- Fqd12/Ftot
QD3f<- Fqd3/Ftot
QD4f<- Fqd4/Ftot
QVf<- Fqv/Ftot

# Resolved in-lever (weighted)
POVw<- POVf*POVdil
PODw<- PODf*PODdil
QD12w<- QD12f*QD12dil
QD3w<- QD3f*QD3dil
QD4w<- QD4f*QD4dil
QVw<- QVf*QVdil

#Resolved jaw adductor in-lever
JAdil<- POVw+PODw+QD12w+QD3w+QD4w+QVw

#Anterior Bite Point Out-lever
ABPdol<- sqrt((JJx-ABPx)^2+(JJy-ABPy)^2+(JJz-ABPz)^2)

#Posterior Bite Point Out-lever
PBPdol<- sqrt((JJx-PBPx)^2+(JJy-PBPy)^2+(JJz-PBPz)^2)

#Anterior Mechanical Advantage (AMA)
AMA<- JAdil/ABPdol

#Posterior Mechanical Advantage (PMA)
PMA<- JAdil/PBPdol

#Anterior Bite Force Production
ABF<- (Ftot*AMA)

#Posterior Bite Force Production
PBF<- (Ftot*PMA)




#Physiological Bite Force Estimates

#Muscle p-CSAs [=(muscle mass*cos theta)/(fiber length*muscle density)]
POVpcsa<- (2.6*1)/(4.483*1.05)
PODpcsa<- (1.2*1)/(4.483*1.05)
QD12pcsa<- (3.5*1)/(4.483*1.05)  #Will likely need to alter
QD3pcsa<- (3.5*1)/(4.483*1.05)
QD4pcsa<- (3.5*1)/(4.483*1.05) 
QVpcsa<- (6.7*1)/(4.483*1.05)
TOTALpcsa<- POVpcsa+PODpcsa+QD12pcsa+QD3pcsa+QD4pcsa+QVpcsa

#Muscle Force Production (Po)--Bilateral (Physiol) 
Fpovp<-POVpcsa*28.9*2
Fpodp<- PODpcsa*28.9*2
Fqd12p<- QD12pcsa*28.9*2
Fqd3p<- QD3pcsa*28.9*2
Fqd4p<- QD4pcsa*14.2*2
Fqvp<- QVpcsa*28.9*2
Ftotp<- Fpovp+Fpodp+Fqd12p+Fqd3p+Fqd4p+Fqvp

#Fraction of total Po
POVfp<- POVpcsa/TOTALpcsa
PODfp<- PODpcsa/TOTALpcsa
QD12fp<- QD12pcsa/TOTALpcsa
QD3fp<- QD3pcsa/TOTALpcsa
QD4fp<- QD4pcsa/TOTALpcsa
QVfp<- QVpcsa/TOTALpcsa

# Resolved in-lever (weighted)
POVwp<- POVfp*POVdil
PODwp<- PODfp*PODdil
QD12wp<- QD12fp*QD12dil
QD3wp<- QD3fp*QD3dil
QD4wp<- QD4fp*QD4dil
QVwp<- QVfp*QVdil

#Resolved jaw adductor in-lever
JAdilp<- POVwp+PODwp+QD12wp+QD3wp+QD4wp+QVwp


#Anterior Bite Force Production
ABFp<- (Ftotp*JAdilp/ABPdol)

#Posterior Bite Force Production
PBFp<- (Ftotp*JAdilp/PBPdol)


