##extraction of sectors##
#define the fraction of the subsidy that is internal
subsidyfrac <-.15
#first the benchmark status quo or complete economy
#let us define the matrix of technical coefficients
#make decimals legible
library(ioanalysis)
#library(OpenMx)
options(scipen = 999) # remove scientific notation
A <- vglad71.IO$A
#record this version of the interaction matrix
Acomp <- A
#and final output (recall that final HH is part of the intermediate sector) 
ftot <- fval
sum(ftot)
sum(ftot - Eval)
I <- diag(x=1,numrows)
L <- solve(I-A)
xcomp <- L%*%ftot
sum(xcomp[1:70])
#check by comparing with the sum of Xval
sum(Xval[1:70])
xcomp[71]
xcomp - A%*%xcomp - ftot
#therefore
#wages
wcomp <- A[71,]*xcomp
xcomp[71]
sum(wcomp)
#surplus
scomp <- A[72,]*xcomp
xcomp[72]
sum(scomp)
#taxes
tcomp <- A[73,]*xcomp
xcomp[73]
sum(tcomp)
xcomp[73]
A[73,21]
#value added
gvacomp <- wcomp + scomp + tcomp 
sum(gvacomp)
sum(xcomp[71:73])
#check shares
wshare <- sum(wcomp)/sum(xcomp[1:70])
wshare
sshare <- sum(scomp)/sum(xcomp[1:70])
sshare
tshare <- sum(tcomp)/sum(xcomp[1:70])
tshare
gvashare <- sum(gvacomp)/sum(xcomp[1:70])
gvashare
wshare/gvashare
gvasharenet <- sum(gvacomp)/(sum(gvacomp)+Totintcomp)
gvasharenet
#exports
ecomp <- Eval

####first extract Al assuming AlOx and El can replace
#contracts with exports####
#scenario Al-full
#total extraction of Al
#reset A
A <- vglad71.IO$A
numrows
A[,21] <- rep(0,numrows)
A[21,] <- rep(0,numrows)
#removal of the subsidy from TaxMinusSubonProd'n
Zval[73,26]
Zval[73,26]+subsidyfrac*174.65
A[73,26]
A[73,26]<- (Zval[73,26]+subsidyfrac*174.65)/Xval[26]
A[73,26]
#record this version of the interaction matrix
AexAl <-A
#reset final demand
ftot <- fval
fval
#now let's set final demand for Al to zero (so everything comes from imports )
ftot[21] <-0
#check on electricity final
ftot[26]
#full replacement of BSL purchases of AlOx with exports
Zval[8,21]
ftot[8] <-ftot[8]+Zval[8,21]
#and full replacement of BSL purchases of
#electricity with exports (inc. sub)
Zval[26,21]
ftot[26] <- ftot[26] + Zval[26,21]+174.65
ftot[26]
#as it is part of imports
#now let's solve the model
L <- solve(I-A)
xexAlfull <- L%*%ftot
#sum up the new output vector
sum(xexAlfull)
sum(xcomp)
#compute the difference
delxexAlfull <- sum(xcomp[1:70] - xexAlfull[1:70])
delxexAlfull
#compute and then display the percentage loss
perchxexAlfull <- delxexAlfull/sum(xcomp[1:70])*100
perchxexAlfull
#wages
wexAlfull <- A[71,]*xexAlfull
xexAlfull[71]
sum(wexAlfull)
#change in wages
#total change in wages
chtotwexAlfull <- xcomp[71]-xexAlfull[71]
#indirect change in wages
chinwexAlfull <- sum(AexAl[71,1:70]*(xcomp[1:70]-xexAlfull[1:70]))
chinwexAlfull
#direct change in wages
chdwexfull <- chtotwexAlfull - chinwexAlfull
chdwexfull
#surplus
sexAlfull <- A[72,]*xexAlfull
xexAlfull[72]
sum(sexAlfull)
#change in surplus
#total change in surplus
chtotsexAlfull <- xcomp[72]-xexAlfull[72]
chtotsexAlfull
#indirect change in surplus
chinsexAlfull <- sum(AexAl[72,1:70]*(xcomp[1:70]-xexAlfull[1:70]))
chinsexAlfull
#direct change in surplus
chdsexfull <- chtotsexAlfull - chinsexAlfull
chdsexfull
#taxes
texAlfull <- A[73,]*xexAlfull
xexAlfull[73]
sum(texAlfull)
A[73,21]
#change in tax
#total change in tax
chtottexAlfull <- xcomp[73]-xexAlfull[73]
#indirect change in tax
chintexAlfull <- sum(AexAl[73,1:70]*(xcomp[1:70]-xexAlfull[1:70]))
chintexAlfull
#direct change in tax
chdtexfull <- chtottexAlfull - chintexAlfull
chdtexfull
#value added
gvaexAlfull <- wexAlfull + sexAlfull + texAlfull 
sum(gvaexAlfull)
sum(xexAlfull[71:73])
# change in value added
sum(gvacomp - gvaexAlfull)
perchgvaexAlfull <- sum(gvacomp - gvaexAlfull)/sum(gvacomp)*100
perchgvaexAlfull
#change in gva
#total change in gva
chtotgvaexAlfull <- chtotwexAlfull + chtotsexAlfull + chtottexAlfull
#indirect change in gva
chingvaexAlfull <- chinwexAlfull + chinsexAlfull + chintexAlfull
chingvaexAlfull
#multiplier
#direct change in gva
chdgvaexAlfull <- chtotgvaexAlfull - chingvaexAlfull
chdgvaexAlfull
#multiplier
chtotgvaexAlfull/chdgvaexAlfull

sum(tcomp - texAlfull)
perchtexAlfull <- sum(tcomp - texAlfull)/sum(tcomp)*100
perchtexAlfull
#indirect impact on wages
A[71,]*(xcomp-xexAlfull)
indwexAlfull <- sum(AexAl[71,1:70]*(xcomp[1:70]-xexAlfull[1:70]))
#total impact on compensation of employees
totwexAlfull <- xcomp[71]-xexAlfull[71]
totwexAlfull
#check direct makes sense
totwexAlfull - indwexAlfull
Zval[71,21]
perchwexAlfull <- (xcomp[71]-xexAlfull[71])/xcomp[71]*100
perchwexAlfull
indperchwexAlfull <- indwexAlfull/xcomp[71]*100
indperchwexAlfull
# changes to exports
eexAlfull <- Eval
eexAlfull[8] <- eexAlfull[8] + Zval[8,21]
eexAlfull[21] <- 0
eexAlfull[26] <- eexAlfull[26] + Zval[26,21]+174.65
#compute the difference
deleexAlfull <- sum(ecomp - eexAlfull)
deleexAlfull
#check
Zval[8,21] - Eval[21] + Zval[26,21]+174.65
#compute and then display the percentage change
sum(eexAlfull[1:70])
sum(ecomp[1:70])
percheexAlfull <- deleexAlfull/sum(ecomp)*100
percheexAlfull
#let's view the relevant part of our new intermediate matrix A using a heatmap
obj = matrix(vglad71.IO$Z,ncol=dim(A)[1])
#plot Zval for the case where Al is extracted
length(xexAlfull)
Id <-diag(1,73)
tr(Id)
#xdiag <- Id[1,]*xexAlfull
#xdiag
#dev.cur()
#obj <- AexAl%*%xdiag
#pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/heatZvalexAl.pdf",width=8, height=6)
#hexAl <- heatmap.io(obj, RS_label, sectors_x = c(5:34,71),
#                    sectors_y = c(5:34,71), FUN = log, max = 7)
#print(plot(hexAl))
#dev.off()


#scenario Al-par-etc
#next the extended extraction cases: total for Al for the other sectors:
# no extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- vglad71.IO$A
#total extraction of Al
A[,21] <- rep(0,numrows)
A[21,] <- rep(0,numrows)
#removal of the subsidy from TaxMinusSubonProd'n
Zval[73,26]
Zval[73,26]+subsidyfrac*174.65
A[73,26]
A[73,26]<- (Zval[73,26]+subsidyfrac*174.65)/Xval[26]
A[73,26]
# 10 percent total extraction of the three potential downstream sectors
# the columns
A[,22] <- .95*A[,22]
A[,23] <- .95*A[,23]
A[,24] <- .95*A[,24]
#the rows
A[22,] <- .95*A[22,]
A[23,] <- .95*A[23,]
A[24,] <- .95*A[24,]
#record this version of the interaction matrix
AexAlpar <-A
#reset f
ftot <- fval
ftot[21] <-0
#full replacement of BSL purchases of AlOx with exports
Zval[8,21]
ftot[8] <-ftot[8]+.95*Zval[8,21]
#and 80% replacement of BSL purchases of electricity
#with exports 
Zval[26,21]
ftot[26] <- ftot[26] + .8*(Zval[26,21]+174.65)
#no need to replace subsidy as it is external here
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.95
ftot[23] <- ftot[23]*.95
ftot[24] <- ftot[24]*.95
#solve the model
L <- solve(I-A)
xexAlpar <- L%*%ftot
#sum up the new output vector
sum(xexAlpar)
#compute the difference
delxexAlpar <- sum(xcomp[1:70] - xexAlpar[1:70])
delxexAlpar
#compute and then display the percentage loss
perchxexAlpar <- delxexAlpar/sum(xcomp[1:70])*100
perchxexAlpar
#this completes the extraction of Al, El and the other sectors
#comparison with Al ie increment over that extraction
perchxexAlpar - perchxexAlfull
#wages
wexAlpar <- A[71,]*xexAlpar
xexAlpar[71]
sum(wexAlpar)
#change in wages
#total change in wages
chtotwexAlpar <- xcomp[71]-xexAlpar[71]
#indirect change in wages
chinwexAlpar <- sum(AexAl[71,1:70]*(xcomp[1:70]-xexAlpar[1:70]))
chinwexAlpar
#direct change in wages
chdwexpar <- chtotwexAlpar - chinwexAlpar
chdwexpar
#surplus
sexAlpar <- A[72,]*xexAlpar
xexAlpar[72]
sum(sexAlpar)
#change in surplus
#total change in surplus
chtotsexAlpar <- xcomp[72]-xexAlpar[72]
chtotsexAlpar
#indirect change in surplus
chinsexAlpar <- sum(AexAl[72,1:70]*(xcomp[1:70]-xexAlpar[1:70]))
chinsexAlpar
#direct change in surplus
chdsexpar <- chtotsexAlpar - chinsexAlpar
chdsexpar
#taxes
texAlpar <- A[73,]*xexAlpar
xexAlpar[73]
sum(texAlpar)
A[73,21]
#change in tax
#total change in tax
chtottexAlpar <- xcomp[73]-xexAlpar[73]
#indirect change in tax
chintexAlpar <- sum(AexAl[73,1:70]*(xcomp[1:70]-xexAlpar[1:70]))
chintexAlpar
#direct change in tax
chdtexpar <- chtottexAlpar - chintexAlpar
chdtexpar
#value added
gvaexAlpar <- wexAlpar + sexAlpar + texAlpar 
sum(gvaexAlpar)
sum(xexAlpar[71:73])
# change in value added
sum(gvacomp - gvaexAlpar)
perchgvaexAlpar <- sum(gvacomp - gvaexAlpar)/sum(gvacomp)*100
perchgvaexAlpar
#change in gva
#total change in gva
chtotgvaexAlpar <- chtotwexAlpar + chtotsexAlpar + chtottexAlpar
#indirect change in gva
chingvaexAlpar <- chinwexAlpar + chinsexAlpar + chintexAlpar
chingvaexAlpar
#multiplier
#direct change in gva
chdgvaexAlpar <- chtotgvaexAlpar - chingvaexAlpar
chdgvaexAlpar
#multiplier
chtotgvaexAlpar/chdgvaexAlpar

#indirect impact on compensation of employees
#A[71,]*(xcomp-xexAlpar)
indwexAlpar <- sum(A[71,]*(xcomp-xexAlpar))
indwexAlpar
#total impact on compensation
totwexAlpar <- xcomp[71]-xexAlpar[71]
totwexAlpar
#check direct makes sense
totwexAlpar - indwexAlpar
#(recall that we have also partially extracted downstream hence 10mm extra)
Zval[71,21]+(1-.95)*(Zval[71,22]+Zval[71,23]+Zval[71,24])
perchwexAlpar <- (xcomp[71]-xexAlpar[71])/xcomp[71]*100
perchwexAlpar
indperchwexAlpar <- indwexAlpar/xcomp[71]*100
indperchwexAlpar
# changes to exports
eexAlpar <- Eval
eexAlpar[8] <- eexAlpar[8] + .95*Zval[8,21]
eexAlpar[21] <- 0
eexAlpar[22] <- eexAlpar[22]*.95
eexAlpar[23] <- eexAlpar[23]*.95
eexAlpar[24] <- eexAlpar[24]*.95
eexAlpar[26] <- eexAlpar[26]+.8*(Zval[26,21]+174.65)
#compute the difference
deleexAlpar <- sum(ecomp - eexAlpar)
deleexAlpar
#check
.95*Zval[8,21] - Eval[21] + .8*(Zval[26,21]+174.65) - (1-.95)*(Eval[22]+Eval[23]+Eval[24])
#compute and then display the percentage change
percheexAlpar <- deleexAlpar/sum(ecomp)*100
percheexAlpar

#scenario Al-none-etc (= Aletc)
#next the extended extraction cases: total for Al for the other sectors:
# no extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- vglad71.IO$A
#total extraction of Al
A[,21] <- rep(0,numrows)
A[21,] <- rep(0,numrows)
#removal of the subsidy from TaxMinusSubonProd'n
Zval[73,26]
Zval[73,26]+subsidyfrac*174.65
A[73,26]
A[73,26]<- (Zval[73,26]+subsidyfrac*174.65)/Xval[26]
A[73,26]
# 10 percent total extraction of the three potential downstream sectors
# the columns
A[,22] <- .95*A[,22]
A[,23] <- .95*A[,23]
A[,24] <- .95*A[,24]
#the rows
A[22,] <- .95*A[22,]
A[23,] <- .95*A[23,]
A[24,] <- .95*A[24,]
#record this version of the interaction matrix
AexAlnone <-A
#reset f
ftot <- fval
ftot[21] <-0
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.95
ftot[23] <- ftot[23]*.95
ftot[24] <- ftot[24]*.95
#solve the model
L <- solve(I-A)
xexAlnone <- L%*%ftot
#sum up the new output vector
sum(xexAlnone)
#compute the difference
delxexAlnone <- sum(xcomp[1:70] - xexAlnone[1:70])
delxexAlnone
#compute and then display the percentage loss
perchxexAlnone <- delxexAlnone/sum(xcomp[1:70])*100
perchxexAlnone
#increment
perchxexAlnone - perchxexAlpar
#wages
wexAlnone <- A[71,]*xexAlnone
xexAlnone[71]
sum(wexAlnone)
#change in wages
#total change in wages
chtotwexAlnone <- xcomp[71]-xexAlnone[71]
#indirect change in wages
chinwexAlnone <- sum(AexAl[71,1:70]*(xcomp[1:70]-xexAlnone[1:70]))
chinwexAlnone
#direct change in wages
chdwexnone <- chtotwexAlnone - chinwexAlnone
chdwexnone
#surplus
sexAlnone <- A[72,]*xexAlnone
xexAlnone[72]
sum(sexAlnone)
#change in surplus
#total change in surplus
chtotsexAlnone <- xcomp[72]-xexAlnone[72]
chtotsexAlnone
#indirect change in surplus
chinsexAlnone <- sum(AexAl[72,1:70]*(xcomp[1:70]-xexAlnone[1:70]))
chinsexAlnone
#direct change in surplus
chdsexnone <- chtotsexAlnone - chinsexAlnone
chdsexnone
#taxes
texAlnone <- A[73,]*xexAlnone
xexAlnone[73]
sum(texAlnone)
A[73,21]
#change in tax
#total change in tax
chtottexAlnone <- xcomp[73]-xexAlnone[73]
#indirect change in tax
chintexAlnone <- sum(AexAl[73,1:70]*(xcomp[1:70]-xexAlnone[1:70]))
chintexAlnone
#direct change in tax
chdtexnone <- chtottexAlnone - chintexAlnone
chdtexnone
#value added
gvaexAlnone <- wexAlnone + sexAlnone + texAlnone 
sum(gvaexAlnone)
sum(xexAlnone[71:73])
# change in value added
sum(gvacomp - gvaexAlnone)
perchgvaexAlnone <- sum(gvacomp - gvaexAlnone)/sum(gvacomp)*100
perchgvaexAlnone
#change in gva
#total change in gva
chtotgvaexAlnone <- chtotwexAlnone + chtotsexAlnone + chtottexAlnone
#indirect change in gva
chingvaexAlnone <- chinwexAlnone + chinsexAlnone + chintexAlnone
chingvaexAlnone
#multiplier
#direct change in gva
chdgvaexAlnone <- chtotgvaexAlnone - chingvaexAlnone
chdgvaexAlnone
#multiplier
chtotgvaexAlnone/chdgvaexAlnone

#indirect impact on compensation of employees
#A[71,]*(xcomp-xexAlnone)
indwexAlnone <- sum(A[71,]*(xcomp-xexAlnone))
indwexAlnone
#total impact on compensation of employees
totwexAlnone <- xcomp[71]-xexAlnone[71]
totwexAlnone
#check direct makes sense
totwexAlnone - indwexAlnone
#(recall that we hve also partially extracted downstream hence 10mm extra)
Zval[71,21]+(1-.95)*(Zval[71,22]+Zval[71,23]+Zval[71,24])
perchwexAlnone <- (xcomp[71]-xexAlnone[71])/xcomp[71]*100
perchwexAlnone
indperchwexAlnone <- indwexAlnone/xcomp[71]*100
indperchwexAlnone
# changes to exports
eexAlnone <- Eval
eexAlnone[21] <- 0
eexAlnone[22] <- eexAlnone[22]*.95
eexAlnone[23] <- eexAlnone[23]*.95
eexAlnone[24] <- eexAlnone[24]*.95
#compute the difference
deleexAlnone <- sum(ecomp - eexAlnone)
deleexAlnone
#check
-Eval[21]- (1-.95)*(Eval[22]+Eval[23]+Eval[24])
#compute and then display the percentage change
percheexAlnone <- deleexAlnone/sum(ecomp)*100
percheexAlnone

##this completes the extraction of Al alone scenarios##
##next Al and El together##

#scenario AlEl-full
#reset A
A <- vglad71.IO$A
#save the hh vue in the interaction matrix
#total extraction of Al
A[,21] <- rep(0,numrows)
A[21,] <- rep(0,numrows)
#removal of the subsidy from TaxMinusSubonProd'n
Zval[73,26]
Zval[73,26]+subsidyfrac*174.65
A[73,26]
A[73,26]<- (Zval[73,26]+subsidyfrac*174.65)/Xval[26]
A[73,26]
#1/3 extraction of El, no need to extract subsidy here
#because it is part of primary inputs (exogenous)
A[,26] <- A[,26]*2/3
A[26,] <- A[26,]*2/3
#adjustment for 192 GPS employees employed at Alsal rate
#ElwnetGPS <- (Zval[71,26]/.117-192)*.117
#ElwnetGPS
#ElwnetGPS/xcomp[26]
#A[71,26]
#A[71,26] <- ElwnetGPS/xcomp[26]
#no coal purchases by GPS
A[6,26] <- 0
#record this interaction matrix
AexAlEl <- A
#reset ftot
ftot <- fval
ftot[21] <-0
#subtract electricity exports and reduce final demand by 1/3
ftot[26]
Eval[26]
ftot[26] <- (ftot[26] -Eval[26])*2/3
ftot[26]
#full replacement of BSL purchases of AlOx with exports
Zval[8,21]
ftot[8] <-ftot[8]+Zval[8,21]
#full replacement of GPS purchases of Coal with exports
Zval[6,26]
ftot[6] <-ftot[6]+Zval[6,26]
#solve the model
L <- solve(I-A)
xexAlElfull <- L%*%ftot
#sum up the new output vector
sum(xexAlElfull)
#compute the difference
delxexAlElfull <- sum(xcomp[1:70] - xexAlElfull[1:70])
delxexAlElfull
#compute and then display the percentage loss
perchxexAlElfull <- delxexAlElfull/sum(xcomp[1:70])*100
perchxexAlElfull
#wages
wexAlElfull <- A[71,]*xexAlElfull
xexAlElfull[71]
sum(wexAlElfull)
#change in wages
#total change in wages
chtotwexAlElfull <- xcomp[71]-xexAlElfull[71]
#indirect change in wages
chinwexAlElfull <- sum(AexAl[71,1:70]*(xcomp[1:70]-xexAlElfull[1:70]))
chinwexAlElfull
#direct change in wages
chdwexfull <- chtotwexAlElfull - chinwexAlElfull
chdwexfull
#surplus
sexAlElfull <- A[72,]*xexAlElfull
xexAlElfull[72]
sum(sexAlElfull)
#change in surplus
#total change in surplus
chtotsexAlElfull <- xcomp[72]-xexAlElfull[72]
chtotsexAlElfull
#indirect change in surplus
chinsexAlElfull <- sum(AexAl[72,1:70]*(xcomp[1:70]-xexAlElfull[1:70]))
chinsexAlElfull
#direct change in surplus
chdsexfull <- chtotsexAlElfull - chinsexAlElfull
chdsexfull
#taxes
texAlElfull <- A[73,]*xexAlElfull
xexAlElfull[73]
sum(texAlElfull)
A[73,21]
#change in tax
#total change in tax
chtottexAlElfull <- xcomp[73]-xexAlElfull[73]
#indirect change in tax
chintexAlElfull <- sum(AexAl[73,1:70]*(xcomp[1:70]-xexAlElfull[1:70]))
chintexAlElfull
#direct change in tax
chdtexfull <- chtottexAlElfull - chintexAlElfull
chdtexfull
#value added
gvaexAlElfull <- wexAlElfull + sexAlElfull + texAlElfull 
sum(gvaexAlElfull)
sum(xexAlElfull[71:73])
# change in value added
sum(gvacomp - gvaexAlElfull)
perchgvaexAlElfull <- sum(gvacomp - gvaexAlElfull)/sum(gvacomp)*100
perchgvaexAlElfull
#change in gva
#total change in gva
chtotgvaexAlElfull <- chtotwexAlElfull + chtotsexAlElfull + chtottexAlElfull
chtotgvaexAlElfull
#indirect change in gva
chingvaexAlElfull <- chinwexAlElfull + chinsexAlElfull + chintexAlElfull
chingvaexAlElfull
#direct change in gva
chdgvaexAlElfull <- chtotgvaexAlElfull - chingvaexAlElfull
chdgvaexAlElfull
#multiplier
multgvaAlElfull <- chtotgvaexAlElfull/chdgvaexAlElfull
multgvaAlElfull
#indirect impact on compensation of employees
#A[71,]*(xcomp-xexAlElfull)
indwexAlElfull <- sum(A[71,]*(xcomp-xexAlElfull))
indwexAlElfull
#total impact on compensation of employees
totwexAlElfull <- xcomp[71]-xexAlElfull[71]
totwexAlElfull
#check direct makes sense
totwexAlElfull - indwexAlElfull
#(recall that we have also partially extracted El)
Zval[71,21]+1/3*Zval[71,26]
perchwexAlElfull <- (xcomp[71]-xexAlElfull[71])/xcomp[71]*100
perchwexAlElfull
indperchwexAlElfull <- indwexAlElfull/xcomp[71]*100
indperchwexAlElfull
#further check on El wages
AexAlEl[71,26]*xexAlElfull[26]/.117
ftot[26]
# changes to exports
eexAlElfull <- Eval
eexAlElfull[6] <- eexAlElfull[6] + Zval[6,26]
eexAlElfull[8] <- eexAlElfull[8] + Zval[8,21]
eexAlElfull[21] <- 0
eexAlElfull[26] <- 0
#compute the difference
deleexAlElfull <- sum(ecomp - eexAlElfull)
deleexAlElfull
#check
Zval[6,26]+ Zval[8,21]-Eval[21]-Eval[26] 
#compute and then display the percentage change
percheexAlElfull <- deleexAlElfull/sum(ecomp)*100
percheexAlElfull

#scenario AlEl-par-etc
#next the extended extraction cases: total for Al and, for the other sectors:
# 70 per cent extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- vglad71.IO$A
#total extraction of Al
A[,21] <- rep(0,numrows)
A[21,] <- rep(0,numrows)
#removal of the subsidy from TaxMinusSubonProd'n
Zval[73,26]
Zval[73,26]+subsidyfrac*174.65
A[73,26]
A[73,26]<- (Zval[73,26]+subsidyfrac*174.65)/Xval[26]
A[73,26]
# 70 percent total extraction of El
A[,26] <- A[,26]*2/3
A[26,] <- A[26,]*2/3
#adjustment for 192 GPS employees employed at Alsal rate
#ElwnetGPS <- (Zval[71,26]/.117-192)*.117
#ElwnetGPS
#ElwnetGPS/xcomp[26]
#A[71,26]
#A[71,26] <- ElwnetGPS/xcomp[26]
#no coal purchases by GPS
A[6,26] <- 0
# 10 percent total extraction of the three potential downstream sectors
# the columns
A[,22] <- .95*A[,22]
A[,23] <- .95*A[,23]
A[,24] <- .95*A[,24]
#the rows
A[22,] <- .95*A[22,]
A[23,] <- .95*A[23,]
A[24,] <- .95*A[24,]
#record this version of the interaction matrix
AexAlElpar <-A
#reset f
ftot <- fval
ftot[21] <-0
#subtract electricity exports and reduce final demand by 1/3
ftot[26] <- (ftot[26] -Eval[26])*2/3
ftot[26]
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.95
ftot[23] <- ftot[23]*.95
ftot[24] <- ftot[24]*.95
#full replacement of BSL purchases of AlOx with exports
Zval[8,21]
ftot[8] <-ftot[8]+.95*Zval[8,21]
#.9 replacement of GPS purchases of Coal with exports
Zval[6,26]
ftot[6] <-ftot[6]+.95*Zval[6,26]
#solve the model
L <- solve(I-A)
xexAlElpar <- L%*%ftot
#sum up the new output vector
sum(xexAlElpar)
#compute the difference
delxexAlElpar <- sum(xcomp[1:70] - xexAlElpar[1:70])
delxexAlElpar
#compute and then display the percentage loss
perchxexAlElpar <- delxexAlElpar/sum(xcomp[1:70])*100
perchxexAlElpar
#wages
wexAlElpar <- A[71,]*xexAlElpar
xexAlElpar[71]
sum(wexAlElpar)
#change in wages
#total change in wages
chtotwexAlElpar <- xcomp[71]-xexAlElpar[71]
#indirect change in wages
chinwexAlElpar <- sum(AexAl[71,1:70]*(xcomp[1:70]-xexAlElpar[1:70]))
chinwexAlElpar
#direct change in wages
chdwexpar <- chtotwexAlElpar - chinwexAlElpar
chdwexpar
#surplus
sexAlElpar <- A[72,]*xexAlElpar
xexAlElpar[72]
sum(sexAlElpar)
#change in surplus
#total change in surplus
chtotsexAlElpar <- xcomp[72]-xexAlElpar[72]
chtotsexAlElpar
#indirect change in surplus
chinsexAlElpar <- sum(AexAl[72,1:70]*(xcomp[1:70]-xexAlElpar[1:70]))
chinsexAlElpar
#direct change in surplus
chdsexpar <- chtotsexAlElpar - chinsexAlElpar
chdsexpar
#taxes
texAlElpar <- A[73,]*xexAlElpar
xexAlElpar[73]
sum(texAlElpar)
A[73,21]
#change in tax
#total change in tax
chtottexAlElpar <- xcomp[73]-xexAlElpar[73]
#indirect change in tax
chintexAlElpar <- sum(AexAl[73,1:70]*(xcomp[1:70]-xexAlElpar[1:70]))
chintexAlElpar
#direct change in tax
chdtexpar <- chtottexAlElpar - chintexAlElpar
chdtexpar
#value added
gvaexAlElpar <- wexAlElpar + sexAlElpar + texAlElpar 
sum(gvaexAlElpar)
sum(xexAlElpar[71:73])
# change in value added
sum(gvacomp - gvaexAlElpar)
perchgvaexAlElpar <- sum(gvacomp - gvaexAlElpar)/sum(gvacomp)*100
perchgvaexAlElpar
#change in gva
#total change in gva
chtotgvaexAlElpar <- chtotwexAlElpar + chtotsexAlElpar + chtottexAlElpar
chtotgvaexAlElpar
#indirect change in gva
chingvaexAlElpar <- chinwexAlElpar + chinsexAlElpar + chintexAlElpar
chingvaexAlElpar
#direct change in gva
chdgvaexAlElpar <- chtotgvaexAlElpar - chingvaexAlElpar
chdgvaexAlElpar
#multiplier
multgvaAlElpar <- chtotgvaexAlElpar/chdgvaexAlElpar
multgvaAlElpar
#indirect impact on compensation of employees
#A[71,]*(xcomp-xexAlElpar)
indwexAlElpar <- sum(A[71,]*(xcomp-xexAlElpar))
indwexAlElpar
#total impact on compensation of employees
totwexAlElpar <- xcomp[71]-xexAlElpar[71]
totwexAlElpar
#check direct makes sense
totwexAlElpar - indwexAlElpar
#(recall that we have also partially extracted downstream)
Zval[71,21]+1/3*Zval[71,26] + (1-.95)*(Zval[71,22]+Zval[71,23]+Zval[71,24])
perchwexAlElpar <- (xcomp[71]-xexAlElpar[71])/xcomp[71]*100
perchwexAlElpar
indperchwexAlElpar <- indwexAlElpar/xcomp[71]*100
indperchwexAlElpar
# changes to exports
eexAlElpar <- Eval
eexAlElpar[6] <- eexAlElpar[6] + .95*Zval[6,26]
eexAlElpar[8] <- eexAlElpar[8] + .95*Zval[8,21]
eexAlElpar[21] <- 0
eexAlElpar[22] <- eexAlElpar[22]*.95
eexAlElpar[23] <- eexAlElpar[23]*.95
eexAlElpar[24] <- eexAlElpar[24]*.95
eexAlElpar[26] <- 0
#compute the difference
deleexAlElpar <- sum(ecomp - eexAlElpar)
deleexAlElpar
#check
.95*Zval[6,26] + .95*Zval[8,21] - Eval[21] - (1-.95)*(Eval[22]+Eval[23]+Eval[24])- Eval[26]
#compute and then display the percentage change
percheexAlElpar <- deleexAlElpar/sum(ecomp)*100
percheexAlElpar

#scenario AlEl-none-etc (=AlElnone)
#next the extended extraction cases: total for Al and, for the other sectors:
# 70 per cent extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- vglad71.IO$A
#total extraction of Al
A[,21] <- rep(0,numrows)
A[21,] <- rep(0,numrows)
#removal of the subsidy from TaxMinusSubonProd'n
Zval[73,26]
Zval[73,26]+subsidyfrac*174.65
A[73,26]
A[73,26]<- (Zval[73,26]+subsidyfrac*174.65)/Xval[26]
A[73,26]
# 70 percent total extraction of El
A[,26] <- A[,26]*2/3
A[26,] <- A[26,]*2/3
#A[26,numrows] <- aEl
#A[26,numrows]
#adjustment for 192 GPS employees employed at Alsal rate
#ElwnetGPS <- (Zval[71,26]/.117-192)*.117
#ElwnetGPS
#ElwnetGPS/xcomp[26]
#A[71,26]
#A[71,26] <- ElwnetGPS/xcomp[26]
#no coal purchases by GPS
A[6,26] <- 0
# 10 percent total extraction of the three potential downstream sectors
# the columns
A[,22] <- .95*A[,22]
A[,23] <- .95*A[,23]
A[,24] <- .95*A[,24]
#the rows
A[22,] <- .95*A[22,]
A[23,] <- .95*A[23,]
A[24,] <- .95*A[24,]
#record this version of the interaction matrix
AexAlElnone <-A
#reset f
ftot <- fval
fval[26]
ftot[26]
ftot[21] <-0
#subtract electricity exports and reduce final demand by 1/3
ftot[26] <- (ftot[26] -Eval[26])*2/3
ftot[26]
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.95
ftot[23] <- ftot[23]*.95
ftot[24] <- ftot[24]*.95
#solve the model
L <- solve(I-A)
xexAlElnone <- L%*%ftot
#sum up the new output vector
sum(xexAlElnone)
#compute the difference
delxexAlElnone <- sum(xcomp[1:70] - xexAlElnone[1:70])
delxexAlElnone
#compute and then display the percentage loss
perchxexAlElnone <- delxexAlElnone/sum(xcomp[1:70])*100
perchxexAlElnone
#wages
wexAlElnone <- A[71,]*xexAlElnone
xexAlElnone[71]
sum(wexAlElnone)
#surplus
sexAlElnone <- A[72,]*xexAlElnone
xexAlElnone[72]
sum(sexAlElnone)
#taxes
texAlElnone <- A[73,]*xexAlElnone
xexAlElnone[73]
sum(texAlElnone)
A[73,21]
#value added
gvaexAlElnone <- wexAlElnone + sexAlElnone + texAlElnone 
sum(gvaexAlElnone)
sum(xexAlElnone[71:73])
# change in value added
sum(gvacomp - gvaexAlElnone)
perchgvaexAlElnone <- sum(gvacomp - gvaexAlElnone)/sum(gvacomp)*100
perchgvaexAlElnone
#indirect impact on compensation of employees
#A[71,]*(xcomp-xexAlElnone)
indwexAlElnone <- sum(A[71,]*(xcomp-xexAlElnone))
indwexAlElnone
#total impact on compensation of employees
totwexAlElnone <- xcomp[71]-xexAlElnone[71]
totwexAlElnone
#check direct makes sense
totwexAlElnone - indwexAlElnone
#vs Zval
Zval[71,21]+ 1/3*Zval[71,26] + (1-.95)*(Zval[71,22]+Zval[71,23]+Zval[71,24])
perchwexAlElnone <- (xcomp[71]-xexAlElnone[71])/xcomp[71]*100
perchwexAlElnone
#quick check
totwexAlElnone/xcomp[71]*100
perchwexAlElnone*xcomp[71]/100
#quick check
totwexAlElnone
indperchwexAlElnone <- indwexAlElnone/xcomp[71]*100
indperchwexAlElnone
# changes to exports
eexAlElnone <- Eval
eexAlElnone[21] <- 0
eexAlElnone[22] <- eexAlElnone[22]*.95
eexAlElnone[23] <- eexAlElnone[23]*.95
eexAlElnone[24] <- eexAlElnone[24]*.95
eexAlElnone[26] <- 0
#compute the difference
deleexAlElnone <- sum(ecomp - eexAlElnone)
deleexAlElnone
#check
- Eval[21] - (1-.95)*(Eval[22]+Eval[23]+Eval[24])- Eval[26]
#compute and then display the percentage change
percheexAlElnone <- deleexAlElnone/sum(ecomp)*100
percheexAlElnone


