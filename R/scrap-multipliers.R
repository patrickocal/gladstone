#truncate to remove Value Added row
vm71io <- vmult71InandOut[-c(71),]
#compare type I and type IIval
vtIlesstII <- vm71io-multInandOut
vtIlesstII[21]
#key sector analysis
#compute the appropriate threshold
keythresh <- median(tIlesstII)
keythresh
keydir2II <- key.sector(vglad71.IO, type =c("direct"),crit = 2 + keythresh)
keydir2II
keytot2II <- key.sector(vglad71.IO, type =c("total"),crit = 2.4)
keytot2II
#
link3 <- linkages(vglad71.IO, regions = 1, type = c("total", "direct"), normalize = FALSE, intra.inter = TRUE)
link3

#Extract the aluminium sector
#first a total extraction
extAltotII <- extraction(vglad71.IO, ES=NULL, regions = 1,sectors = c(21), type="backward.total",aggregate=FALSE,simultaneous=TRUE)
sumAltotII<-sum(extAltotII)
sumAltotII
#next inputs only, so that it is as if the firms keep on importing
extAlinII <- extraction(vglad71.IO, ES=NULL, regions = 1,sectors = c(21), type="backward",aggregate=FALSE,simultaneous=TRUE)
sumAlinII <-sum(extAlinII)
sumAlinII
#difference 
extAltotlessinII <- extAltotII - extAlinII
extAltotlessinII
sum(extAltotlessinII)
#compare with total output
sumx <- sum(X)
sumx
#proportional change
#total
propchextAltotII <- sumAltotII/sumx
propchextAltotII
#only input (backward linkage)
propchextAlinII <- sumAlinII/sumx
propchextAlinII


#a row of ones
ival <- matrix(1,ncol =71)

#Linkages (and also the multipliers in standard vector form)
#output multiplier = Backward linkage, 
#this vector equals the row sum over the columns of the Leontieff matrix L, 
#so for each row of L, we sum the entries. (To get a column, we also transpose.)
BLtval <- t(ival%*%vglad71.IO$L)
BLtbar[21]-BLt[21]
BLtbar[34]-BLt[34]
#input multiplier = Forward linkage,
#this vector equals the column sum over the rows of the Ghoshian matrix G,
#so for each row, we sum the entries.
FLtval <- vglad71.IO$G%*%t(ival)


#"income multiplier (of type II) for sector 1" = $\langle Vbar , L[,1]\rangle$
#vbar <- matrix(as.numeric(unlist(glad70[82, c(4:73,75)])), ncol = dim(Zbar)[1])
#for aluminium
vinc.mult.II.Al <-wbar%*%vglad71.IO$L[,21]/wbar[,21]
vinc.mult.II.Al
vinc.mult.II.Al-inc.mult.I.Al
#for electricity
vinc.mult.II.El <-wbar%*%vglad71.IO$L[,26]/wbar[,26]
vinc.mult.II.El
#for retail
vinc.mult.II.Ret <-wbar%*%vglad71.IO$L[,34]/wbar[,34]
vinc.mult.II.Ret
#and for every sector the whole economy
vinc.mult.II <- wbar%*%vglad71.IO$L/wbar
vinc.mult.II


#"Value-added multiplier (of type II) for sector 21" = $\langle V , L[,21]\rangle$
vval.mult.II.Al <-Vbar%*%vglad71.IO$L[,21]/Vbar[21]
vval.mult.II.Al
vval.mult.II.Al-val.mult.I.Al
#value-added multiplier (of type I) for sector electricity 26 
#and for electricity
vval.mult.II.El <- Vbar%*%vglad71.IO$L[,26]/Vbar[26]
vval.mult.II.El
vval.mult.II.El-val.mult.I.El
