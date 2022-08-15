##extraction of sectors##
#first the benchmark status quo or complete economy
#let us define the matrix of technical coefficients
#make decimals legible
options(scipen = 999) # remove scientific notation
A <- glad71.IO$A
#record this version of the interaction matrix
Acomp <- A
library(ioanalysis)
h <- heatmap.io(Acomp, RS_label, sectors_x = c(5:34,71),
                    sectors_y = c(7:34,71), FUN=log, max = 3)
plot(h)
#and final output (recall that final HH is part of the intermediate sector) 
ftot <- fbar + Ebar
I <- diag(x=1,71)
L <- solve(I-A)
xcomp <- L%*%ftot
sum(xcomp)
#check by comparing with the sum of Xbar
sum(Xbar)
xcomp[71]
xcomp - A%*%xcomp - ftot
#therefore total employee compensation is 
wcomp <- A[71,]*xcomp
sum(A[71,]*xcomp)/sum(xcomp)
#and vue added in the status quo is 
vcomp <- Vbar - wbar + t(Acomp[71,]*xcomp)
sumvcomp <- sum(vcomp)
wbar[71]
vcomp[71]
sum(vcomp)/(sum(xcomp)-sum(Mbar)-xcomp[71]-1000)
sum(xcomp)-sum(Mbar)-xcomp[71]
#compare with original
sum(Vbar)/sum(xcomp)
sum(wbar)/sum(xcomp)
sum(Acomp[71,]*xcomp)/(sum(xcomp)-sum(Mbar))
sum(wbar)/sum(Vbar)
sum(wbar)/12000
sum(wcomp)
ecomp <- Ebar

####first extract Al assuming AlOx and El can replace
#contracts with exports####
#scenario Al-full
#total extraction of Al
#reset A
A <- glad71.IO$A
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#let's view the relevant part of our new intermediate matrix A using a heatmap
obj = matrix(c(A,ftot),ncol=dim(A)[1]+1)
library(ioanalysis)
hexAl <- heatmap.io(obj, RS_label, sectors_x = c(7:34,71),
                    sectors_y = c(7:34,71), FUN=log, max = 7)
plot(hexAl)
#record this version of the interaction matrix
AexAl <-A
#reset final demand
ftot <- fbar + Ebar
#now let's set final demand for Al to zero (so everything comes from imports )
ftot[21] <-0
#check on electricity final
ftot[26]
#full replacement of BSL purchases of AlOx with exports
Zbar[8,21]
ftot[8] <-ftot[8]+Zbar[8,21]
#and full replacement of BSL purchases of electricity with exports
Zbar[26,21]
ftot[26] <- ftot[26] + Zbar[26,21]
#now let's solve the model
L <- solve(I-A)
xexAlfull <- L%*%ftot
#sum up the new output vector
sum(xexAlfull)
#compute the difference
delxexAlfull <- sum(xcomp - xexAlfull)
#compute and then display the percentage loss
perchxexAlfull <- delxexAlfull/sum(xcomp)*100
perchxexAlfull
#indirect impact on compensation of employees
A[71,]*(xcomp-xexAlfull)
indwexAlfull <- sum(A[71,]*(xcomp-xexAlfull))
#total impact on compensation of employees
totwexAlfull <- xcomp[71]-xexAlfull[71]
totwexAlfull
#check direct makes sense
totwexAlfull - indwexAlfull
Zbar[71,21]
perchwexAlfull <- (xcomp[71]-xexAlfull[71])/xcomp[71]*100
perchwexAlfull
indperchwexAlfull <- indwexAlfull/xcomp[71]*100
indperchwexAlfull
#plot Zbar for the case where Al is extracted
xexAlfull
#install.packages("OpenMx")
#library(OpenMx)
xdiag <- vec2diag(xexAlfull)
library(ioanalysis)
h <- heatmap.io(A%*%xdiag, RS_label, sectors_x = c(5:34,71),
                sectors_y = c(5:34,71), FUN=log, max = 7)
plot(h)
#value added
vexAlfull <- Vbar - wbar + t(A[71,]*xexAlfull)
#compute the difference
delvexAlfull <- sum(vcomp - vexAlfull)
#compute and then display the percentage loss
perchvexAlfull <- delvexAlfull/sum(vcomp)*100
perchvexAlfull
# changes to exports
eexAlfull <- Ebar
eexAlfull[8] <- eexAlfull[8] + Zbar[8,21]
eexAlfull[21] <- 0
eexAlfull[26] <- eexAlfull[26] + Zbar[26,21]
#compute the difference
deleexAlfull <- sum(ecomp - eexAlfull)
deleexAlfull
#check
Zbar[8,21] - Ebar[21] + Zbar[26,21]
#compute and then display the percentage change
percheexAlfull <- deleexAlfull/sum(ecomp)*100
percheexAlfull


#scenario Al-par-etc
#next the extended extraction cases: total for Al for the other sectors:
# no extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- glad71.IO$A
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
# 10 percent total extraction of the three potential downstream sectors
A[,22] <- .9*A[,22]
A[,23] <- .9*A[,23]
A[,24] <- .9*A[,24]
#record this version of the interaction matrix
AexAlpar <-A
#reset f
ftot <- fbar + Ebar
ftot[21] <-0
#full replacement of BSL purchases of AlOx with exports
Zbar[8,21]
ftot[8] <-ftot[8]+Zbar[8,21]
#and 80% replacement of BSL purchases of electricity with exports
Zbar[26,21]
ftot[26] <- ftot[26] + .8*Zbar[26,21]
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#solve the model
L <- solve(I-A)
xexAlpar <- L%*%ftot
#sum up the new output vector
sum(xexAlpar)
#compute the difference
delexAlpar <- sum(xcomp - xexAlpar)
#compute and then display the percentage loss
perchxexAlpar <- delexAlpar/sum(xcomp)*100
perchxexAlpar
#this completes the extraction of Al, El and the other sectors
#comparison with Al ie increment over that extraction
perchxexAlpar - perchxexAlfull
#indirect impact on compensation of employees
A[71,]*(xcomp-xexAlpar)
indwexAlpar <- sum(A[71,]*(xcomp-xexAlpar))
indwexAlpar
#total impact on compensation of employees
totwexAlpar <- xcomp[71]-xexAlpar[71]
totwexAlpar
#check direct makes sense
totwexAlpar - indwexAlpar
#(recall that we have also partially extracted downstream hence 10mm extra)
Zbar[71,21]+.1*(Zbar[71,22]+Zbar[71,23]+Zbar[71,24])
perchwexAlpar <- (xcomp[71]-xexAlpar[71])/xcomp[71]*100
perchwexAlpar
indperchwexAlpar <- indwexAlpar/xcomp[71]*100
indperchwexAlpar
#value added
vexAlpar <- Vbar - wbar + t(A[71,]*xexAlpar)
#compute the difference
delvexAlpar <- sum(vcomp - vexAlpar)
#compute and then display the percentage loss
perchvexAlpar <- delvexAlpar/sum(vcomp)*100
perchvexAlpar
# changes to exports
eexAlpar <- Ebar
eexAlpar[8] <- eexAlpar[8] + Zbar[8,21]
eexAlpar[21] <- 0
eexAlpar[22] <- eexAlpar[22]*.9
eexAlpar[23] <- eexAlpar[23]*.9
eexAlpar[24] <- eexAlpar[24]*.9
eexAlpar[26] <- eexAlpar[26]+.8*Zbar[26,21]
#compute the difference
deleexAlpar <- sum(ecomp - eexAlpar)
deleexAlpar
#check
Zbar[8,21] - Ebar[21] + .8*Zbar[26,21] - .1*(Ebar[22]+Ebar[23]+Ebar[24])
#compute and then display the percentage change
percheexAlpar <- deleexAlpar/sum(ecomp)*100
percheexAlpar

#scenario Al-none-etc (= Aletc)
#next the extended extraction cases: total for Al for the other sectors:
# no extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- glad71.IO$A
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
# 10 percent total extraction of the three potential downstream sectors
A[,22] <- .9*A[,22]
A[,23] <- .9*A[,23]
A[,24] <- .9*A[,24]
#record this version of the interaction matrix
AexAlnone <-A
#reset f
ftot <- fbar + Ebar
ftot[21] <-0
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#solve the model
L <- solve(I-A)
xexAlnone <- L%*%ftot
#sum up the new output vector
sum(xexAlnone)
#compute the difference
delexAlnone <- sum(xcomp - xexAlnone)
#compute and then display the percentage loss
perchxexAlnone <- delexAlnone/sum(xcomp)*100
perchxexAlnone
#indirect impact on compensation of employees
A[71,]*(xcomp-xexAlnone)
indwexAlnone <- sum(A[71,]*(xcomp-xexAlnone))
indwexAlnone
#total impact on compensation of employees
totwexAlnone <- xcomp[71]-xexAlnone[71]
totwexAlnone
#check direct makes sense
totwexAlnone - indwexAlnone
#(recall that we have also partially extracted downstream hence 10mm extra)
Zbar[71,21]+.1*(Zbar[71,22]+Zbar[71,23]+Zbar[71,24])
perchwexAlnone <- (xcomp[71]-xexAlnone[71])/xcomp[71]*100
perchwexAlnone
indperchwexAlnone <- indwexAlnone/xcomp[71]*100
indperchwexAlnone
#vue added
vexAlnone <- Vbar - wbar + t(A[71,]*xexAlnone)
#value added
vexAlnone <- Vbar - wbar + t(A[71,]*xexAlnone)
#compute the difference
delvexAlnone <- sum(vcomp - vexAlnone)
#compute and then display the percentage loss
perchvexAlnone <- delvexAlnone/sum(vcomp)*100
perchvexAlnone
# changes to exports
eexAlnone <- Ebar
eexAlnone[21] <- 0
eexAlnone[22] <- eexAlnone[22]*.9
eexAlnone[23] <- eexAlnone[23]*.9
eexAlnone[24] <- eexAlnone[24]*.9
#compute the difference
deleexAlnone <- sum(ecomp - eexAlnone)
deleexAlnone
#check
-Ebar[21]- .1*(Ebar[22]+Ebar[23]+Ebar[24])
#compute and then display the percentage change
percheexAlnone <- deleexAlnone/sum(ecomp)*100
percheexAlnone

##this completes the extraction of Al alone scenarios##
##next Al and El together##

#scenario AlEl-full
#reset A
A <- glad71.IO$A
#save the hh vue in the interaction matrix
aEl <- A[26,71]
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#1/3 extraction of El
A[,26] <- A[,26]*2/3
A[26,] <- A[26,]*2/3
A[26,71] <- aEl
#adjustment for 192 GPS employees employed at Alsal rate
#ElwnetGPS <- (Zbar[71,26]/.117-192)*.117
#ElwnetGPS
#ElwnetGPS/xcomp[26]
#A[71,26]
#A[71,26] <- ElwnetGPS/xcomp[26]
#no coal purchases by GPS
A[6,26] <- 0
#record this interaction matrix
AexAlEl <- A
#reset ftot
ftot <- fbar + Ebar
ftot[21] <-0
#subtract electricity exports and reduce final demand by 1/3
ftot[26] <- (ftot[26] -Ebar[26])*2/3
ftot[26]
#full replacement of BSL purchases of AlOx with exports
Zbar[8,21]
ftot[8] <-ftot[8]+Zbar[8,21]
#full replacement of GPS purchases of Coal with exports
Zbar[6,26]
ftot[6] <-ftot[6]+Zbar[6,26]
#solve the model
L <- solve(I-A)
xexAlElfull <- L%*%ftot
#sum up the new output vector
sum(xexAlElfull)
#compute the difference
delexAlElfull <- sum(xcomp - xexAlElfull)
#compute and then display the percentage loss
perchxexAlElfull <- delexAlElfull/sum(xcomp)*100
perchxexAlElfull
#indirect impact on compensation of employees
A[71,]*(xcomp-xexAlElfull)
indwexAlElfull <- sum(A[71,]*(xcomp-xexAlElfull))
indwexAlElfull
#total impact on compensation of employees
totwexAlElfull <- xcomp[71]-xexAlElfull[71]
totwexAlElfull
#check direct makes sense
totwexAlElfull - indwexAlElfull
#(recall that we have also partially extracted El)
Zbar[71,21]+2/3*Zbar[71,26]
perchwexAlElfull <- (xcomp[71]-xexAlElfull[71])/xcomp[71]*100
perchwexAlElfull
indperchwexAlElfull <- indwexAlElfull/xcomp[71]*100
indperchwexAlElfull
#further check
AexAlEl[71,26]*xexAlElfull[26]
ftot[26]
#vue added
vexAlElfull <- Vbar - wbar + t(A[71,]*xexAlElfull)
#value added
vexAlElfull <- Vbar - wbar + t(A[71,]*xexAlElfull)
#compute the difference
delvexAlElfull <- sum(vcomp - vexAlElfull)
#compute and then display the percentage loss
perchvexAlElfull <- delvexAlElfull/sum(vcomp)*100
perchvexAlElfull
# changes to exports
eexAlElfull <- Ebar
eexAlElfull[6] <- eexAlElfull[6] + Zbar[6,26]
eexAlElfull[8] <- eexAlElfull[8] + Zbar[8,21]
eexAlElfull[21] <- 0
eexAlElfull[26] <- 0
#compute the difference
deleexAlElfull <- sum(ecomp - eexAlElfull)
deleexAlElfull
#check
Zbar[6,26]+ Zbar[8,21]-Ebar[21]-Ebar[26] 
#compute and then display the percentage change
percheexAlElfull <- deleexAlElfull/sum(ecomp)*100
percheexAlElfull

#scenario AlEl-par-etc
#next the extended extraction cases: total for Al and, for the other sectors:
# 70 per cent extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- glad71.IO$A
aEl <- A[26,71]
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
# 70 percent total extraction of El
A[,26] <- A[,26]*2/3
A[26,] <- A[26,]*2/3
A[26,71] <- aEl
#adjustment for 192 GPS employees employed at Alsal rate
#ElwnetGPS <- (Zbar[71,26]/.117-192)*.117
#ElwnetGPS
#ElwnetGPS/xcomp[26]
#A[71,26]
#A[71,26] <- ElwnetGPS/xcomp[26]
#no coal purchases by GPS
A[6,26] <- 0
# 10 percent total extraction of the three potential downstream sectors
A[,22] <- .9*A[,22]
A[,23] <- .9*A[,23]
A[,24] <- .9*A[,24]
#record this version of the interaction matrix
AexAlElpar <-A
#reset f
ftot <- fbar + Ebar
ftot[21] <-0
#subtract electricity exports and reduce final demand by 1/3
ftot[26] <- (ftot[26] -Ebar[26])*2/3
ftot[26]
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#full replacement of BSL purchases of AlOx with exports
Zbar[8,21]
ftot[8] <-ftot[8]+Zbar[8,21]
#.9 replacement of GPS purchases of Coal with exports
Zbar[6,26]
ftot[6] <-ftot[6]+.9*Zbar[6,26]
#solve the model
L <- solve(I-A)
xexAlElpar <- L%*%ftot
#sum up the new output vector
sum(xexAlElpar)
#compute the difference
delexAlElpar <- sum(xcomp - xexAlElpar)
#compute and then display the percentage loss
perchxexAlElpar <- delexAlElpar/sum(xcomp)*100
perchxexAlElpar
#indirect impact on compensation of employees
A[71,]*(xcomp-xexAlElpar)
indwexAlElpar <- sum(A[71,]*(xcomp-xexAlElpar))
indwexAlElpar
#total impact on compensation of employees
totwexAlElpar <- xcomp[71]-xexAlElpar[71]
totwexAlElpar
#check direct makes sense
totwexAlElpar - indwexAlElpar
#(recall that we have also partially extracted downstream)
Zbar[71,21]+2/3*Zbar[71,26] + .1*(Zbar[71,22]+Zbar[71,23]+Zbar[71,24])
perchwexAlElpar <- (xcomp[71]-xexAlElpar[71])/xcomp[71]*100
perchwexAlElpar
indperchwexAlElpar <- indwexAlElpar/xcomp[71]*100
indperchwexAlElpar
#vue added
vexAlElpar <- Vbar - wbar + t(A[71,]*xexAlElpar)
#value added
vexAlElpar <- Vbar - wbar + t(A[71,]*xexAlElpar)
#compute the difference
delvexAlElpar <- sum(vcomp - vexAlElpar)
#compute and then display the percentage loss
perchvexAlElpar <- delvexAlElpar/sum(vcomp)*100
perchvexAlElpar
# changes to exports
eexAlElpar <- Ebar
eexAlElpar[6] <- eexAlElpar[6] + .9*Zbar[6,26]
eexAlElpar[8] <- eexAlElpar[8] + Zbar[8,21]
eexAlElpar[21] <- 0
eexAlElpar[22] <- eexAlElpar[22]*.9
eexAlElpar[23] <- eexAlElpar[23]*.9
eexAlElpar[24] <- eexAlElpar[24]*.9
eexAlElpar[26] <- 0
#compute the difference
deleexAlElpar <- sum(ecomp - eexAlElpar)
deleexAlElpar
#check
.9*Zbar[6,26] + Zbar[8,21] - Ebar[21] - .1*(Ebar[22]+Ebar[23]+Ebar[24])- Ebar[26]
#compute and then display the percentage change
percheexAlElpar <- deleexAlElpar/sum(ecomp)*100
percheexAlElpar

#scenario AlEl-none-etc (=AlElnone)
#next the extended extraction cases: total for Al and, for the other sectors:
# 70 per cent extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- glad71.IO$A
aEl <- A[26,71]
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
# 70 percent total extraction of El
A[,26] <- A[,26]*2/3
A[26,] <- A[26,]*2/3
#A[26,71] <- aEl
#A[26,71]
#adjustment for 192 GPS employees employed at Alsal rate
#ElwnetGPS <- (Zbar[71,26]/.117-192)*.117
#ElwnetGPS
#ElwnetGPS/xcomp[26]
#A[71,26]
#A[71,26] <- ElwnetGPS/xcomp[26]
#no coal purchases by GPS
A[6,26] <- 0
# 10 percent total extraction of the three potential downstream sectors
A[,22] <- .9*A[,22]
A[,23] <- .9*A[,23]
A[,24] <- .9*A[,24]
#record this version of the interaction matrix
AexAlElnone <-A
#reset f
ftot <- fbar + Ebar
fbar[26]
ftot[26]
ftot[21] <-0
#subtract electricity exports and reduce final demand by 1/3
ftot[26] <- (ftot[26] -Ebar[26])*2/3
ftot[26]
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#solve the model
L <- solve(I-A)
xexAlElnone <- L%*%ftot
#sum up the new output vector
sum(xexAlElnone)
#compute the difference
delexAlElnone <- sum(xcomp - xexAlElnone)
#compute and then display the percentage loss
perchxexAlElnone <- delexAlElnone/sum(xcomp)*100
perchxexAlElnone
#indirect impact on compensation of employees
A[71,]*(xcomp-xexAlElnone)
indwexAlElnone <- sum(A[71,]*(xcomp-xexAlElnone))
indwexAlElnone
#total impact on compensation of employees
totwexAlElnone <- xcomp[71]-xexAlElnone[71]
totwexAlElnone
#check direct makes sense
totwexAlElnone - indwexAlElnone
#(recall that we have also partially extracted downstream hence 10mm extra)
#and that we have finely adjusted El labour
#adjustment for 192 GPS employees employed at Alsal rate
Zbar[71,21]+ 2/3*Zbar[71,26] + .1*(Zbar[71,22]+Zbar[71,23]+Zbar[71,24])
perchwexAlElnone <- (xcomp[71]-xexAlElnone[71])/xcomp[71]*100
perchwexAlElnone
#quick check
totwexAlElnone/xcomp[71]*100
perchwexAlElnone*xcomp[71]/100
#quick check
totwexAlElnone
indperchwexAlElnone <- indwexAlElnone/xcomp[71]*100
indperchwexAlElnone
#vue added
vexAlElnone <- Vbar - wbar + t(A[71,]*xexAlElnone)
#value added
vexAlElnone <- Vbar - wbar + t(A[71,]*xexAlElnone)
#compute the difference
delvexAlElnone <- sum(vcomp - vexAlElnone)
#compute and then display the percentage loss
perchvexAlElnone <- delvexAlElnone/sum(vcomp)*100
perchvexAlElnone
# changes to exports
eexAlElnone <- Ebar
eexAlElnone[21] <- 0
eexAlElnone[22] <- eexAlElnone[22]*.9
eexAlElnone[23] <- eexAlElnone[23]*.9
eexAlElnone[24] <- eexAlElnone[24]*.9
eexAlElnone[26] <- 0
#compute the difference
deleexAlElnone <- sum(ecomp - eexAlElnone)
deleexAlElnone
#check
- Ebar[21] - .1*(Ebar[22]+Ebar[23]+Ebar[24])- Ebar[26]
#compute and then display the percentage change
percheexAlElnone <- deleexAlElnone/sum(ecomp)*100
percheexAlElnone

