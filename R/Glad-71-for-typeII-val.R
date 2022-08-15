#Some code to set up an IO table using the R package "ioanalysis"

#first clear the decks
#rm(list=ls())
#then load the data and give it the name "Glad" (when the following commented out, manually)
df <- read.csv("~/OneDrive - The University of Queensland/_QTC/_Air-freight/Modelling/_R/IO-AUS/air-bal-40-Table 1.csv", header=FALSE,row.names=NULL)
#View(df)
df <-AUS114
#the row where the column labels are 
labelrow <- 2
#the column where the row labels are 
labelcol <- 2
#number of intermediate sectors
numint <- 114 
#number of sectors to be included in close out
numrows <-numint+3
numcols <- numint+3
#where to get the data
firstrow <- 4
lastrow <-117
firstcol <- 3
lastcol <- 116


#assuming the table is in the format of the Australia IO table for 2016-2017, the above suffices for the following
#first we set the row numbers for GVA components
wagerow <- 1+lastrow+3 #P1
surrow <-  2+lastrow+3 #P2
taxrow <-  4+lastrow+3 #P4
#then row number for tax on products
tproductrow <- 3+lastrow+3 #P3
#then row number for imports (competing)
mrow <- 6+lastrow+3
#and finally total gva (not currently used)
gvarow <- mrow+4
#similarly for the columns
hhcol <- 1+lastcol+1    #Q1
hhcol
surcol1 <-3+lastcol+1   #Q3
surcol2 <- 6+lastcol+1  #Q6
taxcol1 <- 2+lastcol+1  #Q2
taxcol2 <- 4+lastcol+1  #Q4
taxcol3 <- 5+lastcol+1  #Q5
ecol <- 7+lastcol+1     #Q7
fcol <- 8+lastcol+1     #T5
xcol <- 9+lastcol+1     #T6

#first the components of the IO table
#final consumption
#fEval <- matrix(as.numeric(unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow), c(hhcol,surcol1,surcol2,taxcol1,taxcol2,taxcol3)])), nrow = numrows)
#final hh
fhh <- matrix(as.numeric(gsub(',','',unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow),c(hhcol)]))))
fhh
#final private capital formation
fsur1 <- matrix(as.numeric(gsub(',','',unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow),c(surcol1)]))))
fsur2 <- matrix(as.numeric(gsub(',','',unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow),c(surcol2)]))))
fsur <- fsur1 + fsur2
#final gov 
ftax1 <- matrix(as.numeric(gsub(',','',unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow),c(taxcol1)]))))
ftax2 <- matrix(as.numeric(gsub(',','',unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow),c(taxcol2)]))))
ftax3 <- matrix(as.numeric(gsub(',','',unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow),c(taxcol3)]))))
ftax <- ftax1 + ftax2 + ftax3
#exports
Eval <- matrix(as.numeric(gsub(',','',unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow), ecol]))), nrow = numrows)
Eval
#for ioanalysis
fval <- rep(0,numrows)
#the matrix of intermediate transactions
Zval <- matrix(c(as.numeric(gsub(',','',unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow), c(firstcol:lastcol)])))), numrows, numcols-3)
Zval
Zval <- matrix(c(Zval,fhh,fsur,ftax),numrows,numcols)
Zval
# total supply
Xval <- matrix(as.numeric(gsub(',','',unlist(df[c(firstrow:lastrow,wagerow,surrow,taxrow), xcol]))), nrow = numrows)
Xval
#now for the rows
#value added
vhh <- matrix(as.numeric(gsub(',','',unlist(df[tproductrow,hhcol]))))
vhh
vsur1 <- matrix(as.numeric(gsub(',','',unlist(df[tproductrow,surcol1]))))
vsur1
vsur2 <- matrix(as.numeric(gsub(',','',unlist(df[tproductrow,surcol2]))))
vsur2
vsur = vsur1 + vsur2
vsur
vtax1 <- matrix(as.numeric(gsub(',','',unlist(df[tproductrow,taxcol1]))))
vtax1
vtax2 <- matrix(as.numeric(gsub(',','',unlist(df[tproductrow,taxcol2]))))
vtax2
vtax3 <- matrix(as.numeric(gsub(',','',unlist(df[tproductrow,taxcol3]))))
vtax3
vtax = vtax1 + vtax2 + vtax3
Vval <- matrix(c(as.numeric(gsub(',','',unlist(df[tproductrow, c(firstcol:lastcol)]))),vhh,vsur,vtax), ncol = numcols)
Vval
#Vval <- matrix(as.numeric(unlist(df[c(82), c(firstrow:lastrow,hhcol)])), ncol = 71)
#for the sake of completeness, we also define
mrow
mhh <- matrix(as.numeric(gsub(',','',unlist(df[mrow,hhcol]))))
surcol1
msur1 <- matrix(as.numeric(gsub(',','',unlist(df[mrow,surcol1]))))
msur1
msur2 <- matrix(as.numeric(gsub(',','',unlist(df[mrow,surcol2]))))
msur2
msur = msur1 + msur2
msur
mtax1 <- matrix(as.numeric(gsub(',','',unlist(df[mrow,taxcol1]))))
mtax1
mtax2 <- matrix(as.numeric(gsub(',','',unlist(df[mrow,taxcol2]))))
mtax2
mtax3 <- matrix(as.numeric(gsub(',','',unlist(df[mrow,taxcol3]))))
mtax3
mtax = mtax1 + mtax2 + mtax3
Mval <- matrix(c(as.numeric(gsub(',','',unlist(df[mrow, c(firstcol:lastcol)]))),mhh,msur,mtax), ncol = numcols)
Mval

#indirect taxes (ie GST/VAT on products) and imports for final consumption
#indt <- matrix(as.numeric(df[c(77), c(firstrow:lastrow,hhcol,77,79)]), ncol = 73)
#now for the labels (for sector labels see below)
f_label = df[labelrow,fcol]
E_label = df[labelrow,ecol]
V_label = df[tproductrow,labelcol]
M_label = df[mrow,labelcol]
#fV_label = df[c(77,81),3]
#now the inputoutput object. it seems as if RS_label has to be defined 
#within the object otherwise the dimensions clash (or perhaps this is 
#because it is the only label that needs to have the region defined 
#(hence the two items in the column))
library(ioanalysis)
vglad71.IO <- as.inputoutput(Z=Zval, 
                            RS_label = df[c(firstrow:lastrow,wagerow,surrow,taxrow), c(1,labelcol)],
#                            f=fval, f_label=f_label,
                            E=Eval, E_label=E_label,
                            X=Xval,
                            V=Vval, V_label=V_label,
                            M=Mval, M_label=M_label)
#record value of final: in this case exports
fval<- vglad71.IO$f
fval#fV=fV, fV_label=fV_label)
sum(fval)
sum(vglad71.IO$X[1:numrows])
Totintcomp <- sum(vglad71.IO$A[1:numint,1:numint]%*%vglad71.IO$X[1:numint])
Totintcomp
gva <- sum(vglad71.IO$X[122:124])
gva

mint <- sum(vglad71.IO$M[1:numint])
mint
vint <- sum(vglad71.IO$V[1:numint])
rowkeysum <- Totintcomp + gva + mint + vint
rowkeysum

colkeysum <-sum(Xval[1:numint])
colkeysum
rowkeysum- colkeysum
#check whether the Leontief inverse has negative elements
Lcomp <- vglad71.IO$L
min(Lcomp)
#Locate the worst culprit (if any)
which(Lcomp == min(Lcomp), arr.ind = TRUE)
which(Lcomp <0, arr.ind = TRUE)
Acomp <- vglad71.IO$A
which(Acomp <0, arr.ind = TRUE)
Acomp[73,66]
which(Zval < 0 , arr.ind = TRUE)
Zval[7,72]
Zval[73,21]
sum(Acomp[,21])
sum(Acomp[,26])

dev.off()
#RS_label<-matrix(unlist(as.character(df[firstrow:lastrow, c(1,3)])),nrow=dim(Z)[1])
RS_label = df[c(firstrow:lastrow,wagerow,surrow,taxrow), c(1,3)]
pdf("~/OneDrive - The University of Queensland/_QTC/_Air-freight/Modelling/_R/heatZvalcomp.pdf",width=8, height=6)
h <- heatmap.io(vglad71.IO$Z, RS_label, sectors_x = c(1:numint),
                sectors_y = c(1:numint), FUN=log, max = 7)
dev.off()
print(plot(h))
#input and output multipliers
vmult71InandOut <- multipliers(vglad71.IO,multipliers = "output")
vmult71InandOut
vmult71InandOut[21,]
vmult71InandOut[26]
vmult71InandOut[34]
