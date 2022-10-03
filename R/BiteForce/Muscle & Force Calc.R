#For bull and blacktip sharks
#Distance between two points in 3D

setwd("~/Documents/Shark Project/Digitized Points/")
bite <- read.table("Cleu_Dic01_DigPts.csv", header=TRUE, sep=",")

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

#Perpendicular force of muscle acting on Meckel's cartilage
#using JJ to ABP as reference line for angle of muscle

##POV##
Mmc<- ((JJz-ABPz)/(JJx-ABPx)) #slope of Meckel's cartilage--used as reference
Bmc<- JJz-(Mmc*JJx) #finding intercept for linear equation
Mpov<- ((POViz-POVoz)/(POVix-POVox)) #slope of POV
Bpov<- POViz-(Mpov*POVix)

cm1<- matrix(c(Mpov,-1,Mmc,-1), nrow=2,ncol=2, byrow=TRUE)
rhsmatrix1<- matrix(c(-Bpov,-Bmc), nrow=2,ncol=1, byrow=TRUE) 
inverse1<- solve(cm1)
findint1<- inverse1 %*% rhsmatrix1  #finds point of interesection

Int1x<- findint1[1,1] #creates var for x coord
Int1z<- findint1[2,1] #creates var for z coord

POVbarx<- Int1x-POVox#vector components
POVbarz<- Int1z-POVoz
MCbarx<- Int1x-ABPx
MCbarz<- Int1z-ABPz

POVmag<- sqrt((POVbarx)^2 + (POVbarz)^2)
MCmag<- sqrt((MCbarx)^2 + (MCbarz)^2)

thetaRad1<- acos(((POVbarx*MCbarx)+(POVbarz*MCbarz))/(POVmag*MCmag)) #angle measured in radians
theta1<- thetaRad1*180/pi #angle between Meckel's cartilage and muscle in degrees
phi1<- 90-theta1 #angle away from perpendicular line to MC of given muscle in degrees

Fpovperp<- ((cos(phi1*pi/180))*Fpov)*2 #resultant bilateral force perpendicular to MC

##QD1&2##

Mqd12<- ((QD12iz-QD12oz)/(QD12ix-QD12ox)) #slope of QD12
Bqd12<- QD12iz-(Mqd12*QD12ix)

cm2<- matrix(c(Mqd12,-1,Mmc,-1), nrow=2,ncol=2, byrow=TRUE)
rhsmatrix2<- matrix(c(-Bqd12,-Bmc), nrow=2,ncol=1, byrow=TRUE) 
inverse2<- solve(cm2)
findint2<- inverse2 %*% rhsmatrix2  #finds point of interesection

Int2x<- findint2[1,1] #creates var for x coord
Int2z<- findint2[2,1] #creates var for z coord

QD12barx<- Int2x-QD12ox #vector components
QD12barz<- Int2z-QD12oz
MCbarx<- Int2x-ABPx
MCbarz<- Int2z-ABPz

QD12mag<- sqrt((QD12barx)^2 + (QD12barz)^2)
MCmag<- sqrt((MCbarx)^2 + (MCbarz)^2)

thetaRad2<- acos(((QD12barx*MCbarx)+(QD12barz*MCbarz))/(QD12mag*MCmag))
theta2<- thetaRad2*180/pi #angle between Meckel's cartilage and muscle
phi2<- 90-theta2 #angle away from perpendicular line to MC

Fqd12perp<- ((cos(phi2*pi/180))*Fqd12)*2 #resultant bilateral force perpendicular to MC

##POD##

Mpod<- ((PODiz-PODoz)/(PODix-PODox)) #slope of POD
Bpod<- PODiz-(Mpod*PODix)

cm3<- matrix(c(Mpod,-1,Mmc,-1), nrow=2,ncol=2, byrow=TRUE)
rhsmatrix3<- matrix(c(-Bpod,-Bmc), nrow=2,ncol=1, byrow=TRUE) 
inverse3<- solve(cm3)
findint3<- inverse3 %*% rhsmatrix3  #finds point of interesection

Int3x<- findint3[1,1] #creates var for x coord
Int3z<- findint3[2,1] #creates var for z coord

PODbarx<- Int3x-PODox #vector components
PODbarz<- Int3z-PODoz
MCbarx<- Int3x-ABPx
MCbarz<- Int3z-ABPz

PODmag<- sqrt((PODbarx)^2 + (PODbarz)^2)
MCmag<- sqrt((MCbarx)^2 + (MCbarz)^2)

thetaRad3<- acos(((PODbarx*MCbarx)+(PODbarz*MCbarz))/(PODmag*MCmag))
theta3<- thetaRad3*180/pi #angle between Meckel's cartilage and muscle
phi3<- 90-theta3 #angle away from perpendicular line to MC

Fpodperp<- ((cos(phi3*pi/180))*Fpod)*2 #resultant bilateral force perpendicular to MC

##QD3##

Mqd3<- ((QD3iz-QD3oz)/(QD3ix-QD3ox)) #slope of QD3
Bqd3<- QD3iz-(Mqd3*QD3ix)

cm4<- matrix(c(Mqd3,-1,Mmc,-1), nrow=2,ncol=2, byrow=TRUE)
rhsmatrix4<- matrix(c(-Bqd3,-Bmc), nrow=2,ncol=1, byrow=TRUE) 
inverse4<- solve(cm4)
findint4<- inverse4 %*% rhsmatrix4  #finds point of interesection

Int4x<- findint4[1,1] #creates var for x coord
Int4z<- findint4[2,1] #creates var for z coord

QD3barx<- Int4x-QD3ox #vector components
QD3barz<- Int4z-QD3oz
MCbarx<- Int4x-ABPx
MCbarz<- Int4z-ABPz

QD3mag<- sqrt((QD3barx)^2 + (QD3barz)^2)
MCmag<- sqrt((MCbarx)^2 + (MCbarz)^2)

thetaRad4<- acos(((QD3barx*MCbarx)+(QD3barz*MCbarz))/(QD3mag*MCmag))
theta4<- thetaRad4*180/pi #angle between Meckel's cartilage and muscle
phi4<- 90-theta4 #angle away from perpendicular line to MC

Fqd3perp<- ((cos(phi4*pi/180))*Fqd3)*2 #resultant bilateral force perpendicular to MC

##QD4##

Mqd4<- ((QD4iz-QD4oz)/(QD4ix-QD4ox)) #slope of QD4
Bqd4<- QD4iz-(Mqd4*QD4ix)

cm5<- matrix(c(Mqd4,-1,Mmc,-1), nrow=2,ncol=2, byrow=TRUE)
rhsmatrix5<- matrix(c(-Bqd4,-Bmc), nrow=2,ncol=1, byrow=TRUE) 
inverse5<- solve(cm5)
findint5<- inverse5 %*% rhsmatrix5  #finds point of interesection

Int5x<- findint5[1,1] #creates var for x coord
Int5z<- findint5[2,1] #creates var for z coord

QD4barx<- Int5x-QD4ox #vector components
QD4barz<- Int5z-QD4oz
MCbarx<- Int5x-ABPx
MCbarz<- Int5z-ABPz

QD4mag<- sqrt((QD4barx)^2 + (QD4barz)^2)
MCmag<- sqrt((MCbarx)^2 + (MCbarz)^2)

thetaRad5<- acos(((QD4barx*MCbarx)+(QD4barz*MCbarz))/(QD4mag*MCmag))
theta5<- thetaRad5*180/pi #angle between Meckel's cartilage and muscle
phi5<- 90-theta5 #angle away from perpendicular line to MC

Fqd4perp<- ((cos(phi5*pi/180))*Fqd4)*2 #resultant bilateral force perpendicular to MC

##QV##

Mqv<- ((QViz-QVoz)/(QVix-QVox)) #slope of QV
Bqv<- QViz-(Mqv*QVix)

cm6<- matrix(c(Mqv,-1,Mmc,-1), nrow=2,ncol=2, byrow=TRUE)
rhsmatrix6<- matrix(c(-Bqv,-Bmc), nrow=2,ncol=1, byrow=TRUE) 
inverse6<- solve(cm6)
findint6<- inverse6 %*% rhsmatrix6  #finds point of interesection

Int6x<- findint6[1,1] #creates var for x coord
Int6z<- findint6[2,1] #creates var for z coord

QVbarx<- Int6x-QVox #vector components
QVbarz<- Int6z-QVoz
MCbarx<- Int6x-ABPx
MCbarz<- Int6z-ABPz

QVmag<- sqrt((QVbarx)^2 + (QVbarz)^2)
MCmag<- sqrt((MCbarx)^2 + (MCbarz)^2)

thetaRad6<- acos(((QVbarx*MCbarx)+(QVbarz*MCbarz))/(QVmag*MCmag))
theta6<- thetaRad6*180/pi #angle between Meckel's cartilage and muscle
phi6<- 90-theta6 #angle away from perpendicular line to MC

Fqvperp<- ((cos(phi6*pi/180))*Fqv)*2 #resultant bilateral force perpendicular to MC

#Total force production orthoganol to MC
Ftotperp<- Fpovperp+Fqd12perp+Fpodperp+Fqd3perp+Fqd4perp+Fqvperp

#Anterior Mechanical Advantage (AMA)
AMA<- JAdil/ABPdol

#Posterior Mechanical Advantage (PMA)
PMA<- JAdil/PBPdol

#Anterior Bite Force Production
ABF<- (Ftotperp*AMA)

#Posterior Bite Force Production
PBF<- (Ftotperp*PMA)




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


