
#For bonnetheads and sharks with only 3 QD subdivisions
#Distance between two points in 3D

setwd("/Users/joshcullen/Documents/Shark Project/R Code/")
source("CrossProduct3D.R") #via Chip Hogg in https://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function

setwd("~/Documents/Shark Project/Digitized Points/")
bite <- read.table("BH1_8.1.14_DigPts.csv", header=TRUE, sep=",")

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

QD12ox<- mean(bite[7:9,2]) #QD12 rows may also be just QD1 depending on file
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

QD3ox<- mean(bite[25:27,2]) #QD3 rows may also be QD2 depending on file
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
POVcsa<- bite[1,6] #order of CSAs different in csv files
PODcsa<- bite[19,6]
QD12csa<- bite[7,6]
QD3csa<- bite[25,6]
QD4csa<- bite[31,6]
QVcsa<- bite[13,6]
TOTALcsa<- POVcsa+PODcsa+QD12csa+QD3csa+QD4csa+QVcsa

#Individual Muscle Force Production (Po)--Unilateral 
Fpov<-POVcsa*28.9 #values for red and white muscle spec tension from Lou et al 2002
Fpod<- PODcsa*28.9
Fqd12<- QD12csa*28.9
Fqd3<- QD3csa*28.9
Fqd4<- QD4csa*14.2 
Fqv<- QVcsa*28.9
Ftot<- Fpov+Fpod+Fqd12+Fqd3+Fqd4+Fqv #total unilateral force
Ftotbi<- Ftot*2 # total bilateral force before modifying for angle



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



## Calculation of moments of individual muscles
# Mjj=dotprod(rJJ*crossprod(pv*fv))
# BF= 2*sum(Mjj)/out lever

#dummy axis of jaw joints
JJopp<- c(JJx, (JJy-15), JJz) #right JJ y-coord (-15) chosen at random
rJJ<- (JJopp-JJ)/15 #unit vector of JJ axis

#muscle unit vectors (insertion coord-origin coord)/magnitude
rPOV<- (POVi-POVo)/POVdm
rQD12<- (QD12i-QD12o)/QD12dm
rQD3<- (QD3i-QD3o)/QD3dm
rQV<- (QVi-QVo)/QVdm
rPOD<- (PODi-PODo)/PODdm
rQD4<- (QD4i-QD4o)/QD4dm  ###Double-check that component vectors make sense###

#creation of muscle force vectors (unit vector*force)
fvPOV<- rPOV*Fpov
fvQD12<- rQD12*Fqd12
fvQD3<- rQD3*Fqd3
fvQV<- rQV*Fqv
fvPOD<- rPOD*Fpod
fvQD4<- rQD4*Fqd4

#calculating position vector of muscle insertion from JJ
pvPOV<- POVi-JJ
pvQD12<- QD12i-JJ
pvQD3<- QD3i-JJ
pvQV<- QVi-JJ
pvPOD<- PODi-JJ
pvQD4<- QD4i-JJ

#cross prod of force vector and position vector followed by dot prod with rJJ
Mpov<- sum(rJJ*CrossProduct3D(fvPOV,pvPOV))
Mqd12<- sum(rJJ*CrossProduct3D(fvQD12,pvQD12))
Mqd3<- sum(rJJ*CrossProduct3D(fvQD3,pvQD3))
Mqv<- sum(rJJ*CrossProduct3D(fvQV,pvQV))  
Mpod<- sum(rJJ*CrossProduct3D(fvPOD,pvPOD))
Mqd4<- sum(rJJ*CrossProduct3D(fvQD4,pvQD4))

Mtot<- sum(Mpov+Mqd12+Mqd3+Mqv+Mpod+Mqd4)

#create perpendicular vector to that of ABP/PBP-JJ-JJopp
rABP<- t(as.matrix(ABP-JJ)) #vector from JJ to ABP
rPBP<- t(as.matrix(PBP-JJ)) #vector from JJ to PBP
uvABP<-rABP/sqrt(rABP[,1]^2+rABP[,2]^2+rABP[,3]^2) #unit vector
uvPBP<-rPBP/sqrt(rPBP[,1]^2+rPBP[,2]^2+rPBP[,3]^2) #unit vector

rABF<- t(as.matrix(CrossProduct3D(uvABP,rJJ))) #cross product of unit vectors results in perp vector to describe BF direction
rPBF<- t(as.matrix(CrossProduct3D(uvPBP,rJJ))) # for PBF direction

#BF production
ABFnew<- Mtot/rABF[1,] #summed moments divded by unit vector perp to JJ-BP plane
PBFnew<- Mtot/rPBF[1,]

###last calc for BF is incorrect and is only unilateral; how to incorporate perp vector to BP-JJ plane?
