#Some code to set up an IO table using the R package "ioanalysis"

#first clear the decks
rm(list=ls())
#then load the data and give it the name "Glad" (when the following commented out, manually)
#Glad <- read.csv(
#  "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_R/IO-70-canonical/Glad70-numbers.csv",
#  header=FALSE)
#then shorten the name and have two viewable versions of the data (useful for large datasets)
glad70 <- Glad
library(ioanalysis)
#View(glad70)
wagerow <-75
hhcol <- 75
#first the components of the IO table
#starting with the matrix of intermediate transactions
Zbar <- matrix(as.numeric(unlist(glad70[c(4:73,wagerow), c(4:73,hhcol)])), ncol = 71)
#then the columns
#final consumption
fEhhbar <- matrix(as.numeric(unlist(glad70[c(4:73,wagerow), 82])), nrow = dim(Zbar)[1])
#In this case, final consumption includes exports
Ebar <- matrix(as.numeric(unlist(glad70[c(4:73,wagerow), 81])), nrow = dim(Zbar)[1])
# and final hh consumption,
hhbar <- matrix(as.numeric(unlist(glad70[c(4:73,wagerow), hhcol])), nrow = dim(Zbar)[1])
#so exports need to be subtracted to get net final consumption 
fbar <- fEhhbar-Ebar-hhbar
# total supply
Xbar <- matrix(as.numeric(unlist(glad70[c(4:73,wagerow), 83])), nrow = dim(Zbar)[1])
#now for the rows
#value added
Vbar <- matrix(as.numeric(unlist(glad70[c(82), c(4:73,hhcol)])), ncol = dim(Zbar)[1])
#value added includes employee compensation "wbar"
wbar <- matrix(as.numeric(unlist(glad70[c(75), c(4:73,hhcol)])), ncol = dim(Zbar)[1])
#imports
Mbar <- matrix(as.numeric(unlist(glad70[79, c(4:73,hhcol)])), ncol = dim(Zbar)[1])
#indirect taxes (ie GST/VAT on products) and imports for final consumption
#fV <- matrix(as.numeric(glad70[c(77,79), c(82)]), nrow = 2)
#now for the labels (for sector labels see below)
f_label = glad70[3,82]
E_label = glad70[3,81]
V_label = glad70[82,3]
M_label = glad70[79,3]
#fV_label = glad70[c(77,81),3]
#now the inputoutput object. it seems as if RS_label has to be defined 
#within the object otherwise the dimensions clash (or perhaps this is 
#because it is the only label that needs to have the region defined 
#(hence the two items in the column))
library(ioanalysis)
glad71.IO <- as.inputoutput(Z=Zbar, 
                          RS_label = glad70[c(4:73,75), c(1,3)],
                          f=fbar, f_label=f_label,
                          E=Ebar, E_label=E_label,
                          X=Xbar,
                          V=Vbar, V_label=V_label,
                          M=Mbar, M_label=M_label)
                          #fV=fV, fV_label=fV_label)
#RS_label<-matrix(unlist(as.character(glad70[4:73, c(1,3)])),nrow=dim(Z)[1])
RS_label = glad70[c(4:73,wagerow,surrow,taxrow), c(1,3)]
h <- heatmap.io(glad71.IO$Z, RS_label, sectors_x = c(5:34,71),
                sectors_y = c(5:34,71), FUN=log, max = 7)
plot(h)
pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/heatZcomp.pdf",height=6,width=8)
#input and output multipliers
mult71InandOut <-multipliers(glad71.IO,multipliers = c("input","output"))
mult71InandOut[21]
mult71InandOut[26]
mult71InandOut[34]
#truncate to remove Value Added row
m71io <- mult71InandOut[-c(71),]
#compare type I and type II
tIlesstII <- m71io-multInandOut
tIlesstII[21]
#key sector analysis
#compute the appropriate threshold
keythresh <- median(tIlesstII)
keythresh
keydir2II <- key.sector(glad71.IO, type =c("direct"),crit = 2 + keythresh)
keydir2II
keytot2II <- key.sector(glad71.IO, type =c("total"),crit = 2.4)
keytot2II
#
link3 <- linkages(glad71.IO, regions = 1, type = c("total", "direct"), normalize = FALSE, intra.inter = TRUE)
link3

#Extract the aluminium sector
#first a total extraction
extAltotII <- extraction(glad71.IO, ES=NULL, regions = 1,sectors = c(21), type="backward.total",aggregate=FALSE,simultaneous=TRUE)
sumAltotII<-sum(extAltotII)
sumAltotII
#next inputs only, so that it is as if the firms keep on importing
extAlinII <- extraction(glad71.IO, ES=NULL, regions = 1,sectors = c(21), type="backward",aggregate=FALSE,simultaneous=TRUE)
sumAlinII <-sum(extAlinII)
sumAlinII
#difference 
extAltotlessinII <- extAltotII - extAlinII
extAltotlessinII
sum(extAltotlessinII)
#compare with total output
#proportional change
#total
propchextAltotII <- sumAltotII/sumx
propchextAltotII
#only input (backward linkage)
propchextAlinII <- sumAlinII/sumx
propchextAlinII


#a row of ones
ibar <- matrix(1,ncol =71)

#Linkages (and also the multipliers in standard vector form)
#output multiplier = Backward linkage, 
#this vector equals the row sum over the columns of the Leontieff matrix L, 
#so for each row of L, we sum the entries. (To get a column, we also transpose.)
BLtbar <- t(ibar%*%glad71.IO$L)
BLtbar
BLtbar[21]-BLt[21]
BLtbar[34]-BLt[34]
#input multiplier = Forward linkage,
#this vector equals the column sum over the rows of the Ghoshian matrix G,
#so for each row, we sum the entries.
FLtbar <- glad71.IO$G%*%t(ibar)
FLtbar

#"income multiplier (of type II) for sector 1" = $\langle Vbar , L[,1]\rangle$
#vbar <- matrix(as.numeric(unlist(glad70[82, c(4:73,75)])), ncol = dim(Zbar)[1])
#for aluminium
inc.mult.II.Al <-wbar%*%glad71.IO$L[,21]/wbar[,21]
inc.mult.II.Al
inc.mult.II.Al-inc.mult.I.Al
#for electricity
inc.mult.II.El <-wbar%*%glad71.IO$L[,26]/wbar[,26]
inc.mult.II.El
#for retail
inc.mult.II.Ret <-wbar%*%glad71.IO$L[,34]/wbar[,34]
inc.mult.II.Ret
#and for every sector the whole economy
inc.mult.II <- wbar%*%glad71.IO$L/wbar
inc.mult.II


#"Value-added multiplier (of type II) for sector 21" = $\langle V , L[,21]\rangle$
val.mult.II.Al <-Vbar%*%glad71.IO$L[,21]/Vbar[21]
val.mult.II.Al
val.mult.II.Al-val.mult.I.Al
#value-added multiplier (of type I) for sector electricity 26 
#and for electricity
val.mult.II.El <- Vbar%*%glad71.IO$L[,26]/Vbar[26]
val.mult.II.El
val.mult.II.El-val.mult.I.El

