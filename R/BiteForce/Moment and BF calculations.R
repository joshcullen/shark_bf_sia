## Calculation of moments of individual muscles
# Mjj=dotprod(rJJ*crossprod(pv*fv))
# BF= 2*sum(Mjj)/out lever

#dummy axis of jaw joints
JJopp<- c(JJx, (JJy-15), JJz) #right JJ y coord chosen at random
rJJ<- JJopp-JJ #unit vector of J axis

#muscle unit vectors (insertion coord-origin coord)/magnitude
rPOV<- (POVi-POVo)/POVdm
rQD12<- (QD12i-QD12o)/QD12dm
rQV<- (QVi-QVo)/QVdm
rPOD<- (PODi-PODo)/PODdm
rQD3<- (QD3i-QD3o)/QD3dm
rQD4<- (QD4i-QD4o)/QD4dm  ###Double-check that component vectors make sense###

#creation of muscle force vectors (unit vector*force)
fvPOV<- rPOV*Fpov
fvQD12<- rQD12*Fqd12
fvQV<- rQV*Fqv
fvPOD<- rPOD*Fpod
fvQD3<- rQD3*Fqd3
fvQD4<- rQD4*Fqd4

#calculating position vector of muscle insertion from JJ
pvPOV<- JJ-POVi
pvQD12<- JJ-QD12i
pvQV<- JJ-QVi
pvPOD<- JJ-PODi
pvQD3<- JJ-QD3i
pvQD4<- JJ-QD4i

#cross prod of force vector and position vector followed by dot prod with rJJ
Mpov<- sum(rJJ*(fvPOV %*% pvPOV))
Mqd12<- sum(rJJ*(fvQD12 %*% pvQD12))
Mqv<- sum(rJJ*(fvQV %*% pvQV))
Mpod<- sum(rJJ*(fvPOD %*% pvPOD))
Mqd3<- sum(rJJ*(fvQD3 %*% pvQD3))
Mqd4<- sum(rJJ*(fvQD4 %*% pvQD4))

Mtot<- sum(Mpov+Mqd12+Mqv+Mpod+Mqd3+Mqd4)

#create perpendicular unit vector to that of ABP/PBP-JJ-JJopp
rABP<- ABP-JJ #unit vector from JJ to ABP
rPBP<- PBP-JJ #unit vector from JJ to PBP

rABF<- rABP %*% rJJ #cross product of unit vectors results in perp unit vector to describe BF direction
rPBF<- rPBP %*% rJJ # for PBF direction

#BF production
ABFnew<- Mtot/rABF[1,] #summed moments divded by unit vector perp to JJ-BP plane
PBFnew<- Mtot/rPBF[1,]