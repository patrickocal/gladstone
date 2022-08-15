##extraction of sectors##
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
ftot
sum(ftot)
sum(ftot - Eval)
I <- diag(x=1,numrows)
L <- solve(I-A)
xcomp <- L%*%ftot
sum(xcomp[1:numint])
#check by comparing with the sum of Xval
sum(Xval[1:numint])
xcomp - A%*%xcomp - ftot
#therefore
#wages
wcomp <- A[eval(numint+1),]*xcomp
xcomp[eval(numint+1)]
sum(wcomp)
#surplus
scomp <- A[eval(numint+2),]*xcomp
xcomp[eval(numint+2)]
sum(scomp)
#taxes
tcomp <- A[eval(numint+3),]*xcomp
xcomp[eval(numint+3)]
sum(tcomp)
xcomp[eval(numint+3)]
#value added
gvacomp <- wcomp + scomp + tcomp 
sum(gvacomp)
sum(xcomp[eval(numint+1):eval(numint+3)])
#check shares
wshare <- sum(wcomp)/sum(xcomp[1:numint])
wshare
sshare <- sum(scomp)/sum(xcomp[1:numint])
sshare
tshare <- sum(tcomp)/sum(xcomp[1:numint])
tshare
gvashare <- sum(gvacomp)/sum(xcomp[1:numint])
gvashare
wshare/gvashare
gvasharenet <- sum(gvacomp)/(sum(gvacomp)+Totintcomp)
gvasharenet
#exports
ecomp <- ftot
sum(ecomp)
####first extract Al assuming AlOx and El can replace
#contracts with exports####
#scenario Al-full
#total extraction of Al
#reset A
A <- vglad71.IO$A
numrows
A[,1:7] <- A[,1:7]*.8
A[1:7,8:numrows] <- A[1:7,8:numrows]*.8
A
#record this version of the interaction matrix
AexAl <-A
#reset final demand
ftot <- fval
fval
#now let's set final demand for Al to zero (so everything comes from imports )
ftot[1:7] <- ftot[1:7]*.8
ftot
#now let's solve the model
L <- solve(I-A)
xexAlfull <- L%*%ftot
#sum up the new output vector
sum(xexAlfull)
sum(xcomp)
#compute the difference
delxexAlfull <- sum(xcomp[1:numint] - xexAlfull[1:numint])
delxexAlfull
#compute and then display the percentage loss
perchxexAlfull <- delxexAlfull/sum(xcomp[1:numint])*100
perchxexAlfull
#wages
wexAlfull <- A[eval(numint+1),]*xexAlfull
xexAlfull[eval(numint+1)]
sum(wexAlfull)
#change in wages
#total change in wages
chtotwexAlfull <- xcomp[eval(numint+1)]-xexAlfull[eval(numint+1)]
chtotwexAlfull
#indirect change in wages
chinwexAlfull <- sum(AexAl[eval(numint+1),1:numint]*(xcomp[1:numint]-xexAlfull[1:numint]))
chinwexAlfull
#direct change in wages
chdwexAlfull <- chtotwexAlfull - chinwexAlfull
chdwexAlfull
#surplus
sexAlfull <- A[eval(numint+2),]*xexAlfull
xexAlfull[eval(numint+2)]
sum(sexAlfull)
#change in surplus
#total change in surplus
chtotsexAlfull <- xcomp[eval(numint+2)]-xexAlfull[eval(numint+2)]
chtotsexAlfull
#indirect change in surplus
chinsexAlfull <- sum(AexAl[eval(numint+2),1:numint]*(xcomp[1:numint]-xexAlfull[1:numint]))
chinsexAlfull
#direct change in surplus
chdsexAlfull <- chtotsexAlfull - chinsexAlfull
chdsexAlfull
#taxes
texAlfull <- A[eval(numint+3),]*xexAlfull
xexAlfull[eval(numint+3)]
sum(texAlfull)
#change in tax
#total change in tax
chtottexAlfull <- xcomp[eval(numint+3)]-xexAlfull[eval(numint+3)]
#indirect change in tax
chintexAlfull <- sum(AexAl[eval(numint+3),1:numint]*(xcomp[1:numint]-xexAlfull[1:numint]))
chintexAlfull
#direct change in tax
chdtexAlfull <- chtottexAlfull - chintexAlfull
chdtexAlfull
#value added
gvaexAlfull <- wexAlfull + sexAlfull + texAlfull 
sum(gvaexAlfull)
sum(xexAlfull[eval(numint+1):eval(numint+3)])
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
#check
chdwexAlfull + chdsexAlfull + chdtexAlfull
#multiplier
multgvaAlfull <- chtotgvaexAlfull/chdgvaexAlfull
multgvaAlfull
#total impact on compensation of employees
totwexAlfull <- xcomp[eval(numint+1)]-xexAlfull[eval(numint+1)]
totwexAlfull
#check direct makes sense
totwexAlfull - chinwexAlfull
Zval[eval(numint+1),21]
perchwexAlfull <- (xcomp[eval(numint+1)]-xexAlfull[eval(numint+1)])/xcomp[eval(numint+1)]*100
perchwexAlfull
indperchwexAlfull <- chinwexAlfull/xcomp[eval(numint+1)]*100
indperchwexAlfull
# changes to exports
eexAlfull <- ftot
eexAlfull
#compute the difference
deleexAlfull <- sum(ecomp - eexAlfull)
deleexAlfull
#check
sum(ecomp[1:7])*(1-.8)
#compute and then display the percentage change
sum(eexAlfull[1:numint])
sum(ecomp[1:numint])
percheexAlfull <- deleexAlfull/sum(ecomp)*100
percheexAlfull
#let's view the relevant part of our new intermediate matrix A using a heatmap
#obj = matrix(vglad71.IO$Z,ncol=dim(A)[1])
#heatmap(obj)
#plot Zval for the case where Al is extracted
#length(xexAlfull)

#xdiag <- Id[1,]*xexAlfull
#xdiag
#dev.cur()
#obj <- AexAl%*%xdiag
#pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/heatZvalexAl.pdf",width=8, height=6)
#hexAl <- heatmap.io(obj, RS_label, sectors_x = c(5:34,eval(numint+1)),
#                    sectors_y = c(5:34,eval(numint+1)), FUN = log, max = 7)
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
#and 80% replacement of BSL purchases of electricity with exports
Zval[26,21]
ftot[26] <- ftot[26] + .8*Zval[26,21]
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
delxexAlpar <- sum(xcomp[1:numint] - xexAlpar[1:numint])
delxexAlpar
#compute and then display the percentage loss
perchxexAlpar <- delxexAlpar/sum(xcomp[1:numint])*100
perchxexAlpar
#this completes the extraction of Al, El and the other sectors
#comparison with Al ie increment over that extraction
perchxexAlpar - perchxexAlfull
#wages
wexAlpar <- A[eval(numint+1),]*xexAlpar
xexAlpar[eval(numint+1)]
sum(wexAlpar)
#change in wages
#total change in wages
chtotwexAlpar <- xcomp[eval(numint+1)]-xexAlpar[eval(numint+1)]
#indirect change in wages
chinwexAlpar <- sum(AexAlpar[eval(numint+1),1:numint]*(xcomp[1:numint]-xexAlpar[1:numint]))
chinwexAlpar
#direct change in wages
chdwexAlpar <- chtotwexAlpar - chinwexAlpar
chdwexAlpar
#surplus
sexAlpar <- A[eval(numint+2),]*xexAlpar
xexAlpar[eval(numint+2)]
sum(sexAlpar)
#change in surplus
#total change in surplus
chtotsexAlpar <- xcomp[eval(numint+2)]-xexAlpar[eval(numint+2)]
chtotsexAlpar
#indirect change in surplus
chinsexAlpar <- sum(AexAlpar[eval(numint+2),1:numint]*(xcomp[1:numint]-xexAlpar[1:numint]))
chinsexAlpar
#direct change in surplus
chdsexAlpar <- chtotsexAlpar - chinsexAlpar
chdsexAlpar
#taxes
texAlpar <- A[eval(numint+3),]*xexAlpar
xexAlpar[eval(numint+3)]
sum(texAlpar)
A[eval(numint+3),21]
#change in tax
#total change in tax
chtottexAlpar <- xcomp[eval(numint+3)]-xexAlpar[eval(numint+3)]
#indirect change in tax
chintexAlpar <- sum(AexAlpar[eval(numint+3),1:numint]*(xcomp[1:numint]-xexAlpar[1:numint]))
chintexAlpar
#direct change in tax
chdtexAlpar <- chtottexAlpar - chintexAlpar
chdtexAlpar
#value added
gvaexAlpar <- wexAlpar + sexAlpar + texAlpar 
sum(gvaexAlpar)
sum(xexAlpar[eval(numint+1):eval(numint+3)])
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
#check
chdwexAlpar + chdsexAlpar + chdtexAlpar
#multiplier
multgvaAlpar <- chtotgvaexAlpar/chdgvaexAlpar
multgvaAlpar
#indirect impact on compensation of employees
#A[eval(numint+1),]*(xcomp-xexAlpar)
indwexAlpar <- sum(A[eval(numint+1),]*(xcomp-xexAlpar))
indwexAlpar
#total impact on compensation
totwexAlpar <- xcomp[eval(numint+1)]-xexAlpar[eval(numint+1)]
totwexAlpar
#check direct makes sense
totwexAlpar - indwexAlpar
#(recall that we have also partially extracted downstream hence 10mm extra)
Zval[eval(numint+1),21]+(1-.95)*(Zval[eval(numint+1),22]+Zval[eval(numint+1),23]+Zval[eval(numint+1),24])
perchwexAlpar <- (xcomp[eval(numint+1)]-xexAlpar[eval(numint+1)])/xcomp[eval(numint+1)]*100
perchwexAlpar
indperchwexAlpar <- indwexAlpar/xcomp[eval(numint+1)]*100
indperchwexAlpar
# changes to exports
eexAlpar <- Eval
eexAlpar[8] <- eexAlpar[8] + .95*Zval[8,21]
eexAlpar[21] <- 0
eexAlpar[22] <- eexAlpar[22]*.95
eexAlpar[23] <- eexAlpar[23]*.95
eexAlpar[24] <- eexAlpar[24]*.95
eexAlpar[26] <- eexAlpar[26]+.8*Zval[26,21]
#compute the difference
deleexAlpar <- sum(ecomp - eexAlpar)
deleexAlpar
#check
.95*Zval[8,21] - Eval[21] + .8*Zval[26,21] - (1-.95)*(Eval[22]+Eval[23]+Eval[24])
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
delxexAlnone <- sum(xcomp[1:numint] - xexAlnone[1:numint])
delxexAlnone
#compute and then display the percentage loss
perchxexAlnone <- delxexAlnone/sum(xcomp[1:numint])*100
perchxexAlnone
#increment
perchxexAlnone - perchxexAlpar
#wages
wexAlnone <- A[eval(numint+1),]*xexAlnone
xexAlnone[eval(numint+1)]
sum(wexAlnone)
#change in wages
#total change in wages
chtotwexAlnone <- xcomp[eval(numint+1)]-xexAlnone[eval(numint+1)]
#indirect change in wages
chinwexAlnone <- sum(AexAlnone[eval(numint+1),1:numint]*(xcomp[1:numint]-xexAlnone[1:numint]))
chinwexAlnone
#direct change in wages
chdwexAlnone <- chtotwexAlnone - chinwexAlnone
chdwexAlnone
#surplus
sexAlnone <- A[eval(numint+2),]*xexAlnone
xexAlnone[eval(numint+2)]
sum(sexAlnone)
#change in surplus
#total change in surplus
chtotsexAlnone <- xcomp[eval(numint+2)]-xexAlnone[eval(numint+2)]
chtotsexAlnone
#indirect change in surplus
chinsexAlnone <- sum(AexAlnone[eval(numint+2),1:numint]*(xcomp[1:numint]-xexAlnone[1:numint]))
chinsexAlnone
#direct change in surplus
chdsexAlnone <- chtotsexAlnone - chinsexAlnone
chdsexAlnone
#taxes
texAlnone <- A[eval(numint+3),]*xexAlnone
xexAlnone[eval(numint+3)]
sum(texAlnone)
A[eval(numint+3),21]
#change in tax
#total change in tax
chtottexAlnone <- xcomp[eval(numint+3)]-xexAlnone[eval(numint+3)]
#indirect change in tax
chintexAlnone <- sum(AexAlnone[eval(numint+3),1:numint]*(xcomp[1:numint]-xexAlnone[1:numint]))
chintexAlnone
#direct change in tax
chdtexAlnone <- chtottexAlnone - chintexAlnone
chdtexAlnone
#value added
gvaexAlnone <- wexAlnone + sexAlnone + texAlnone 
sum(gvaexAlnone)
sum(xexAlnone[eval(numint+1):eval(numint+3)])
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
#check
chdwexAlnone + chdsexAlnone + chdtexAlnone
#multiplier
multgvaAlnone <- chtotgvaexAlnone/chdgvaexAlnone
multgvaAlnone
#indirect impact on compensation of employees
#A[eval(numint+1),]*(xcomp-xexAlnone)
indwexAlnone <- sum(A[eval(numint+1),]*(xcomp-xexAlnone))
indwexAlnone
#total impact on compensation of employees
totwexAlnone <- xcomp[eval(numint+1)]-xexAlnone[eval(numint+1)]
totwexAlnone
#check direct makes sense
totwexAlnone - indwexAlnone
#(recall that we hve also partially extracted downstream hence 10mm extra)
Zval[eval(numint+1),21]+(1-.95)*(Zval[eval(numint+1),22]+Zval[eval(numint+1),23]+Zval[eval(numint+1),24])
perchwexAlnone <- (xcomp[eval(numint+1)]-xexAlnone[eval(numint+1)])/xcomp[eval(numint+1)]*100
perchwexAlnone
indperchwexAlnone <- indwexAlnone/xcomp[eval(numint+1)]*100
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
#1/3 extraction of El
A[,26] <- A[,26]*2/3
A[26,] <- A[26,]*2/3
#adjustment for 192 GPS employees employed at Alsal rate
#ElwnetGPS <- (Zval[eval(numint+1),26]/.117-192)*.117
#ElwnetGPS
#ElwnetGPS/xcomp[26]
#A[eval(numint+1),26]
#A[eval(numint+1),26] <- ElwnetGPS/xcomp[26]
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
delxexAlElfull <- sum(xcomp[1:numint] - xexAlElfull[1:numint])
delxexAlElfull
#compute and then display the percentage loss
perchxexAlElfull <- delxexAlElfull/sum(xcomp[1:numint])*100
perchxexAlElfull
#wages
wexAlElfull <- A[eval(numint+1),]*xexAlElfull
xexAlElfull[eval(numint+1)]
sum(wexAlElfull)
#change in wages
#total change in wages
chtotwexAlElfull <- xcomp[eval(numint+1)]-xexAlElfull[eval(numint+1)]
#indirect change in wages
chinwexAlElfull <- sum(AexAlEl[eval(numint+1),1:numint]*(xcomp[1:numint]-xexAlElfull[1:numint]))
chinwexAlElfull
#direct change in wages
chdwexpar <- chtotwexAlElfull - chinwexAlElfull
chdwexpar
#surplus
sexAlElfull <- A[eval(numint+2),]*xexAlElfull
xexAlElfull[eval(numint+2)]
sum(sexAlElfull)
#change in surplus
#total change in surplus
chtotsexAlElfull <- xcomp[eval(numint+2)]-xexAlElfull[eval(numint+2)]
chtotsexAlElfull
#indirect change in surplus
chinsexAlElfull <- sum(AexAlEl[eval(numint+2),1:numint]*(xcomp[1:numint]-xexAlElfull[1:numint]))
chinsexAlElfull
#direct change in surplus
chdsexAlElfull <- chtotsexAlElfull - chinsexAlElfull
chdsexAlElfull
#taxes
texAlElfull <- A[eval(numint+3),]*xexAlElfull
xexAlElfull[eval(numint+3)]
sum(texAlElfull)
A[eval(numint+3),21]
#change in tax
#total change in tax
chtottexAlElfull <- xcomp[eval(numint+3)]-xexAlElfull[eval(numint+3)]
#indirect change in tax
chintexAlElfull <- sum(AexAlEl[eval(numint+3),1:numint]*(xcomp[1:numint]-xexAlElfull[1:numint]))
chintexAlElfull
#direct change in tax
chdtexpar <- chtottexAlElfull - chintexAlElfull
chdtexpar
#value added
gvaexAlElfull <- wexAlElfull + sexAlElfull + texAlElfull 
sum(gvaexAlElfull)
sum(xexAlElfull[eval(numint+1):eval(numint+3)])
# change in value added
sum(gvacomp - gvaexAlElfull)
perchgvaexAlElfull <- sum(gvacomp - gvaexAlElfull)/sum(gvacomp)*100
perchgvaexAlElfull
#change in gva
#total change in gva
chtotgvaexAlElfull <- chtotwexAlElfull + chtotsexAlElfull + chtottexAlElfull
#indirect change in gva
chingvaexAlElfull <- chinwexAlElfull + chinsexAlElfull + chintexAlElfull
chingvaexAlElfull
#multiplier
#direct change in gva
chdgvaexAlElfull <- chtotgvaexAlElfull - chingvaexAlElfull
chdgvaexAlElfull
#multiplier
multgvaAlElfull <- chtotgvaexAlElfull/chdgvaexAlElfull
#indirect impact on compensation of employees
#A[eval(numint+1),]*(xcomp-xexAlElfull)
indwexAlElfull <- sum(A[eval(numint+1),]*(xcomp-xexAlElfull))
indwexAlElfull
#total impact on compensation of employees
totwexAlElfull <- xcomp[eval(numint+1)]-xexAlElfull[eval(numint+1)]
totwexAlElfull
#check direct makes sense
totwexAlElfull - indwexAlElfull
#(recall that we have also partially extracted El)
Zval[eval(numint+1),21]+1/3*Zval[eval(numint+1),26]
perchwexAlElfull <- (xcomp[eval(numint+1)]-xexAlElfull[eval(numint+1)])/xcomp[eval(numint+1)]*100
perchwexAlElfull
indperchwexAlElfull <- indwexAlElfull/xcomp[eval(numint+1)]*100
indperchwexAlElfull
#further check on El wages
AexAlEl[eval(numint+1),26]*xexAlElfull[26]/.117
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
# numint per cent extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- vglad71.IO$A
#total extraction of Al
A[,21] <- rep(0,numrows)
A[21,] <- rep(0,numrows)
# numint percent total extraction of El
A[,26] <- A[,26]*2/3
A[26,] <- A[26,]*2/3
#adjustment for 192 GPS employees employed at Alsal rate
#ElwnetGPS <- (Zval[eval(numint+1),26]/.117-192)*.117
#ElwnetGPS
#ElwnetGPS/xcomp[26]
#A[eval(numint+1),26]
#A[eval(numint+1),26] <- ElwnetGPS/xcomp[26]
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
delxexAlElpar <- sum(xcomp[1:numint] - xexAlElpar[1:numint])
delxexAlElpar
#compute and then display the percentage loss
perchxexAlElpar <- delxexAlElpar/sum(xcomp[1:numint])*100
perchxexAlElpar
#wages
wexAlElpar <- A[eval(numint+1),]*xexAlElpar
xexAlElpar[eval(numint+1)]
sum(wexAlElpar)
#change in wages
#total change in wages
chtotwexAlElpar <- xcomp[eval(numint+1)]-xexAlElpar[eval(numint+1)]
#indirect change in wages
chinwexAlElpar <- sum(AexAlElpar[eval(numint+1),1:numint]*(xcomp[1:numint]-xexAlElpar[1:numint]))
chinwexAlElpar
#direct change in wages
chdwexpar <- chtotwexAlElpar - chinwexAlElpar
chdwexpar
#surplus
sexAlElpar <- A[eval(numint+2),]*xexAlElpar
xexAlElpar[eval(numint+2)]
sum(sexAlElpar)
#change in surplus
#total change in surplus
chtotsexAlElpar <- xcomp[eval(numint+2)]-xexAlElpar[eval(numint+2)]
chtotsexAlElpar
#indirect change in surplus
chinsexAlElpar <- sum(AexAlElpar[eval(numint+2),1:numint]*(xcomp[1:numint]-xexAlElpar[1:numint]))
chinsexAlElpar
#direct change in surplus
chdsexAlElpar <- chtotsexAlElpar - chinsexAlElpar
chdsexAlElpar
#taxes
texAlElpar <- A[eval(numint+3),]*xexAlElpar
xexAlElpar[eval(numint+3)]
sum(texAlElpar)
A[eval(numint+3),21]
#change in tax
#total change in tax
chtottexAlElpar <- xcomp[eval(numint+3)]-xexAlElpar[eval(numint+3)]
#indirect change in tax
chintexAlElpar <- sum(AexAlElpar[eval(numint+3),1:numint]*(xcomp[1:numint]-xexAlElpar[1:numint]))
chintexAlElpar
#direct change in tax
chdtexpar <- chtottexAlElpar - chintexAlElpar
chdtexpar
#value added
gvaexAlElpar <- wexAlElpar + sexAlElpar + texAlElpar 
sum(gvaexAlElpar)
sum(xexAlElpar[eval(numint+1):eval(numint+3)])
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
#A[eval(numint+1),]*(xcomp-xexAlElpar)
indwexAlElpar <- sum(A[eval(numint+1),]*(xcomp-xexAlElpar))
indwexAlElpar
#total impact on compensation of employees
totwexAlElpar <- xcomp[eval(numint+1)]-xexAlElpar[eval(numint+1)]
totwexAlElpar
#check direct makes sense
totwexAlElpar - indwexAlElpar
#(recall that we have also partially extracted downstream)
Zval[eval(numint+1),21]+1/3*Zval[eval(numint+1),26] + (1-.95)*(Zval[eval(numint+1),22]+Zval[eval(numint+1),23]+Zval[eval(numint+1),24])
perchwexAlElpar <- (xcomp[eval(numint+1)]-xexAlElpar[eval(numint+1)])/xcomp[eval(numint+1)]*100
perchwexAlElpar
indperchwexAlElpar <- indwexAlElpar/xcomp[eval(numint+1)]*100
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
# numint per cent extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- vglad71.IO$A
#total extraction of Al
A[,21] <- rep(0,numrows)
A[21,] <- rep(0,numrows)
# numint percent total extraction of El
A[,26] <- A[,26]*2/3
A[26,] <- A[26,]*2/3
#A[26,numrows] <- aEl
#A[26,numrows]
#adjustment for 192 GPS employees employed at Alsal rate
#ElwnetGPS <- (Zval[eval(numint+1),26]/.117-192)*.117
#ElwnetGPS
#ElwnetGPS/xcomp[26]
#A[eval(numint+1),26]
#A[eval(numint+1),26] <- ElwnetGPS/xcomp[26]
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
delxexAlElnone <- sum(xcomp[1:numint] - xexAlElnone[1:numint])
delxexAlElnone
#compute and then display the percentage loss
perchxexAlElnone <- delxexAlElnone/sum(xcomp[1:numint])*100
perchxexAlElnone
#wages
wexAlElnone <- A[eval(numint+1),]*xexAlElnone
xexAlElnone[eval(numint+1)]
sum(wexAlElnone)
#change in wages
#total change in wages
chtotwexAlElnone <- xcomp[eval(numint+1)]-xexAlElnone[eval(numint+1)]
#indirect change in wages
chinwexAlElnone <- sum(AexAlElnone[eval(numint+1),1:numint]*(xcomp[1:numint]-xexAlElnone[1:numint]))
chinwexAlElnone
#direct change in wages
chdwexpar <- chtotwexAlElnone - chinwexAlElnone
chdwexpar
#surplus
sexAlElnone <- A[eval(numint+2),]*xexAlElnone
xexAlElnone[eval(numint+2)]
sum(sexAlElnone)
#change in surplus
#total change in surplus
chtotsexAlElnone <- xcomp[eval(numint+2)]-xexAlElnone[eval(numint+2)]
chtotsexAlElnone
#indirect change in surplus
chinsexAlElnone <- sum(AexAlElnone[eval(numint+2),1:numint]*(xcomp[1:numint]-xexAlElnone[1:numint]))
chinsexAlElnone
#direct change in surplus
chdsexAlElnone <- chtotsexAlElnone - chinsexAlElnone
chdsexAlElnone
#taxes
texAlElnone <- A[eval(numint+3),]*xexAlElnone
xexAlElnone[eval(numint+3)]
sum(texAlElnone)
A[eval(numint+3),21]
#change in tax
#total change in tax
chtottexAlElnone <- xcomp[eval(numint+3)]-xexAlElnone[eval(numint+3)]
#indirect change in tax
chintexAlElnone <- sum(AexAlElnone[eval(numint+3),1:numint]*(xcomp[1:numint]-xexAlElnone[1:numint]))
chintexAlElnone
#direct change in tax
chdtexpar <- chtottexAlElnone - chintexAlElnone
chdtexpar
#value added
gvaexAlElnone <- wexAlElnone + sexAlElnone + texAlElnone 
sum(gvaexAlElnone)
sum(xexAlElnone[eval(numint+1):eval(numint+3)])
# change in value added
sum(gvacomp - gvaexAlElnone)
perchgvaexAlElnone <- sum(gvacomp - gvaexAlElnone)/sum(gvacomp)*100
perchgvaexAlElnone
#change in gva
#total change in gva
chtotgvaexAlElnone <- chtotwexAlElnone + chtotsexAlElnone + chtottexAlElnone
#indirect change in gva
chingvaexAlElnone <- chinwexAlElnone + chinsexAlElnone + chintexAlElnone
chingvaexAlElnone
#multiplier
#direct change in gva
chdgvaexAlElnone <- chtotgvaexAlElnone - chingvaexAlElnone
chdgvaexAlElnone
#multiplier
multgvaAlElnone <- chtotgvaexAlElnone/chdgvaexAlElnone
#indirect impact on compensation of employees
#A[eval(numint+1),]*(xcomp-xexAlElnone)
indwexAlElnone <- sum(A[eval(numint+1),]*(xcomp-xexAlElnone))
indwexAlElnone
#total impact on compensation of employees
totwexAlElnone <- xcomp[eval(numint+1)]-xexAlElnone[eval(numint+1)]
totwexAlElnone
#check direct makes sense
totwexAlElnone - indwexAlElnone
#vs Zval
Zval[eval(numint+1),21]+ 1/3*Zval[eval(numint+1),26] + (1-.95)*(Zval[eval(numint+1),22]+Zval[eval(numint+1),23]+Zval[eval(numint+1),24])
perchwexAlElnone <- (xcomp[eval(numint+1)]-xexAlElnone[eval(numint+1)])/xcomp[eval(numint+1)]*100
perchwexAlElnone
#quick check
totwexAlElnone/xcomp[eval(numint+1)]*100
perchwexAlElnone*xcomp[eval(numint+1)]/100
#quick check
totwexAlElnone
indperchwexAlElnone <- indwexAlElnone/xcomp[eval(numint+1)]*100
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


