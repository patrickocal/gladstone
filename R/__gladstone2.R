#Some code to set up an IO table using the R package "ioanalysis"
#first clear the decks
rm(list=ls())
#options(scipen = 999) # remove scientific notation
#then load the data (when this is commented out, manually)
#Glad70.numbers <- read.csv(
#  "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_R/_Glad70-all/_Glad-complete-mm.csv",
#  header=FALSE)
#import tailored sector labels for larger diagrams
#RS_label.dat <- read.csv("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_R/_Glad70-all/70sector-labels.csv")
#then shorten the name and have two viewable versions of the data (useful for large datasets)
glad70 <- Glad

#first the components of the IO table
#starting with the matrix of intermediate transactions
Z <- matrix(as.numeric(unlist(glad70[4:116, 4:116])), ncol = 70)
#then the columns
#final consumption
f <- matrix(as.numeric(unlist(glad70[4:73, c(82)])), nrow = dim(Z)[1])
#exports
E <- matrix(as.numeric(unlist(glad70[4:73, 81])), nrow = dim(Z)[1])
#if f includes exports, the latter needs to be subtracted to get net final consumption 
f <- f-E
# total supply
X <- matrix(as.numeric(unlist(glad70[4:73, 83])), nrow = dim(Z)[1])
#now for the rows
#value added
V <- matrix(as.numeric(unlist(glad70[c(82), 4:73])), ncol = dim(Z)[1])
#wages
w <- matrix(as.numeric(unlist(glad70[c(75), 4:73])), ncol = dim(Z)[1])
#imports
M <- matrix(as.numeric(unlist(glad70[79, 4:73])), ncol = dim(Z)[1])
#the matrix of final demand's value added
fV <- matrix(as.numeric(unlist(glad70[c(82), c(82)])), nrow = 1)
#now for the labels (for sector labels see below)
f_label <- matrix(as.character(unlist(glad70[3,82])), nrow = 1)
E_label <- matrix(as.character(unlist(glad70[3,81])), nrow = 1)
V_label <- matrix(as.character(unlist(glad70[c(82),3])), nrow = 1)
M_label <- matrix(as.character(unlist(glad70[79,3])), nrow =1)
fV_label <- matrix(as.character(unlist(glad70[82,3])), nrow = 1)
#now the inputoutput object. it seems as if RS_label has to be defined 
#within the object otherwise the dimensions clash (or perhaps this is 
#because it is the only label that needs to have the region defined 
#(hence the two items in the column))
glad.IO <- as.inputoutput(Z=Z, 
                          RS_label = glad70[4:73, c(1,3)],
                          f=f, f_label=f_label,
                          E=E, E_label=E_label,
                          X=X,
                          V=V, V_label=V_label,
                          M=M, M_label=M_label,
                          fV=fV, fV_label=fV_label)
RS_label<-matrix(as.character(unlist(glad70[4:73, c(1,3)])),nrow=dim(Z)[1])
#RS_label = RS_label.dat[4:73, c(1,3)]
hZ <- heatmap.io(glad.IO$Z, RS_label, FUN=log, max = 8)
plot(hZ)
hL <- heatmap.io(glad.IO$L, RS_label, FUN=log, max = 3)
plot(hL)
hG <- heatmap.io(glad.IO$G, RS_label, FUN=log, max = 3)
plot(hG)
hB <- heatmap.io(glad.IO$B, RS_label, FUN=log, max = 3)
plot(hB)
#type I input and output multipliers
multInandOut <-multipliers(glad.IO,multipliers = c("input","output"))
multInandOut[21]
multInandOut[26]
multInandOut[34]
#income multiplier
M1 <- multipliers(glad.IO, multipliers = "wage", wage.row = 1)
M1
#key sector analysis
keydir2 <- key.sector(glad.IO, type =c("direct"),crit = 2)
keydir2
keytot2 <- key.sector(glad.IO, type =c("total"),crit = 2)
keytot2
#
link3 <- linkages(glad.IO, regions = 1, type = c("total", "direct"), normalize = FALSE, intra.inter = TRUE)
link3
#Extract the aluminium sector
extAlonly = extraction(glad.IO, ES=NULL, regions = 1,sectors = c(21), type="backward.total",aggregate=TRUE,simultaneous=TRUE)
extAlonly
write.table(extAlonly, "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_R/_Glad70-all/extAlonly.csv", sep=",")
#Extract just electricity
extElonly = extraction(glad.IO, ES=NULL, regions = 1,sectors = c(26), type="backward.total",aggregate=TRUE,simultaneous=TRUE)
extElonly
write.table(extElonly, "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_R/_Glad70-all/extElonly.csv", sep=",")
#Extract both aluminium and electricity simulaneously
extAlsimEl = extraction(glad.IO, ES=NULL, regions = 1,sectors = c(21,26), type="backward.total",aggregate=TRUE,simultaneous=TRUE)
extAlsimEl
write.table(extAlsimEl, "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_R/_Glad70-all/extAlsimEl.csv", sep=",")
#Extract aluminium and electricity sequentially
extAlthenEl = extraction(glad.IO, ES=NULL, regions = 1,sectors = c(21,26), type="backward.total",aggregate=TRUE,simultaneous=FALSE)
extAlthenEl
write.table(extAlthenEl, "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_R/_Glad70-all/extAlthenEl.csv", sep=",")

#a row of ones
i <- matrix(1,ncol =70)

#Linkages (and also the multipliers in standard vector form)
#output multiplier = Backward linkage, 
#this vector equals the row sum over the columns of the Leontieff matrix L, 
#so for each row of L, we sum the entries. (To get a column, we also transpose.)
BLt <- t(i%*%glad.IO$L)
t(BLt)
BLt[34]
BLt[21]
#input multiplier = Forward linkage,
#this vector equals the column sum over the rows of the Ghoshian matrix G,
#so for each row, we sum the entries.
FLt <- glad.IO$G%*%t(i)
t(FLt)
#"income multiplier (of type I) for sector Al" = $\langle V , L[,21]\rangle$
inc.mult.I.Al <-w%*%glad.IO$L[,21]/w[21]
inc.mult.I.Al
w[21]/sum(w)
#"income multiplier (of type I) for sector electricity " = $\langle w , L[,21]\rangle$
inc.mult.I.El <-w%*%glad.IO$L[,26]/w[26]
inc.mult.I.El
#for retail
inc.mult.I.Ret <-w%*%glad.IO$L[,34]/w[,34]
inc.mult.I.Ret


#and for every sector the whole economy
inc.mult.I <- w%*%glad.IO$L/w
inc.mult.I

#"Value-added multiplier (of type I) for sector 21" = $\langle V , L[,21]\rangle$
val.mult.I.Al <-V%*%glad.IO$L[,21]/V[21]
val.mult.I.Al
V[21]/sum(V)
#value-added multiplier (of type I) for sector electricity 26 
#and for electricity
val.mult.I.El <- V%*%glad.IO$L[,26]/V[26]
val.mult.I.El
V[26]/sum(V)
#glad.IO$L%*%glad.IO$f
#and for every sector the whole economy
val.mult.I <- V%*%glad.IO$L/V
val.mult.I
glad.IO$L[,21]




