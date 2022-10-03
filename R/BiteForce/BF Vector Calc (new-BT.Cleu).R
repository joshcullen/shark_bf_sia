#Need to use CrossProduct3D fxn to calculate properly

setwd("/Users/joshcullen/Documents/Shark Project/R Code/")
source("CrossProduct3D.R") #via Chip Hogg in https://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function

rABF<- CrossProduct3D(rABP,rJJ)
as.matrix(rABF) #allows for indexing
rABF<- t(rABF) #transpose from col to row vector
uvABF<- rABF/sqrt(rABF[1,1]^2+rABF[1,2]^2+rABF[1,3]^2) #creates unit vector

rPBF<-CrossProduct3D(rPBP,rJJ)
as.matrix(rPBF)
rPBF<- t(rPBF)
uvPBF<-rPBF/sqrt(rPBF[,1]^2+rPBF[,2]^2+rPBF[,3]^2)


#Proj(fv)=(fv dot uvABF)*uvABF for all muscle fvs

POVabfv<- sum(fvPOV*uvABF)*uvABF #coords for new force vector
PODabfv<- sum(fvPOD*uvABF)*uvABF
QD1abfv<- sum(fvQD1*uvABF)*uvABF
QD2abfv<- sum(fvQD2*uvABF)*uvABF
QD3abfv<- sum(fvQD3*uvABF)*uvABF
QD4abfv<- sum(fvQD4*uvABF)*uvABF
QVabfv<- sum(fvQV*uvABF)*uvABF

POVpbfv<- sum(fvPOV*uvPBF)*uvPBF #coords for new force vector
PODpbfv<- sum(fvPOD*uvPBF)*uvPBF
QD1pbfv<- sum(fvQD1*uvPBF)*uvPBF
QD2pbfv<- sum(fvQD2*uvPBF)*uvPBF
QD3pbfv<- sum(fvQD3*uvPBF)*uvPBF
QD4pbfv<- sum(fvQD4*uvPBF)*uvPBF
QVpbfv<- sum(fvQV*uvPBF)*uvPBF

POVabf<- sum(fvPOV*uvABF) #transformed muscle forces orthog to JJ plane
PODabf<- sum(fvPOD*uvABF)
QD1abf<- sum(fvQD1*uvABF)
QD2abf<- sum(fvQD2*uvABF)
QD3abf<- sum(fvQD3*uvABF)
QD4abf<- sum(fvQD4*uvABF)
QVabf<- sum(fvQV*uvABF)

POVpbf<- sum(fvPOV*uvPBF)
PODpbf<- sum(fvPOD*uvPBF)
QD1pbf<- sum(fvQD1*uvPBF)
QD2pbf<- sum(fvQD2*uvPBF)
QD3pbf<- sum(fvQD3*uvPBF)
QD4pbf<- sum(fvQD4*uvPBF)
QVpbf<- sum(fvQV*uvPBF)

ABF2<- (POVabf+PODabf+QD1abf+QD2abf+QD3abf+QD4abf+QVabf)*JAdil/ABPdol*2 #bilateral ABF
PBF2<- (POVpbf+PODpbf+QD1pbf+QD2pbf+QD3pbf+QD4pbf+QVpbf)*JAdil/PBPdol*2 #bilateral PBF
