#then load the data and give it the name "Glad" (when the following commented out, manually)
table8 <- data.table::data.table(readr::read_csv("sio4-19-preRAS.csv"))
table8 <- readxl::read_excel("520905500108.xlsx",
                             sheet = "Table 8",
                             col_names=TRUE)
table8 <- data.table::data.table(table8)
#some basic tidying up and label assignment
#colnames(sio4) <- as.character(sio4[1,])
#sio4 <- sio4[-1,]
#first change the name of ourcode71
#colnames(sio4)[which(colnames(sio4)=="ourcode71")] <- "our71names"
#in preparation for ioanalysis

#sio4 <- dplyr::mutate(sio4, state = "qld")
#sio4$state <- rep("",length(sio4$our71names))
#sio4$state <- dplyr::if_else(is.na(sio4$state),
#                      "qld",
#                      "qld")
View(table8)
The location of things in number terms.
```{r message=TRUE, warning=TRUE, include=FALSE}
#the row where the column labels are 
#number of intermediate sectors
numint <- 71 
#number of sectors to be closed out
numclose <- 1
#number of sectors to be included in close out
numrows <-numint+numclose
numcols <- numint+numclose
#where to get the data



firstrow <- which(sio4$our71names == "01 Cattle & calves")
lastrow <-which(sio4$our71names == "71 Other services")
firstcol <- which(colnames(sio4) == "01 Cattle & calves")
lastcol <- which(colnames(sio4) == "71 Other services")
ucol <- which(colnames(sio4) == "U")
sio4$U
#first the final demand columns
hhcol <- which(colnames(sio4)== "H")    #Q1
r1col <- which(colnames(sio4)== "R1")
r2col <- which(colnames(sio4)== "R2")
ecol <- which(colnames(sio4)== "E")     #Q7
#fcol <-      #T5
xcol <- which(colnames(sio4) == "X")     #T6

#assuming the table is in the format of the Australia IO table for 2016-2017, the above suffices for the following
#first we set the row numbers for GVA components
wagerow <- which(sio4$our71names == "W") #P1
wagerow
surrow <-  which(sio4$our71names == "S") #P2
taxrow <-  which(sio4$our71names == "Tv") #P4
#then row number for tax on products
tcrow <- which(sio4$our71names == "Tc") #P3
tcrow
#then row number for imports (competing)
mrow <- which(sio4$our71names == "M")
mrow
vrow <- which(sio4$our71names== "V")
vrow

rrow <- which(sio4$our71names == "R")
rrow
#and finally total gva (not currently used)
#gvarow <- mrow+4

```
 The final demand components of the IO table
```{r}
#first
#final consumption
#fEval <- matrix(as.numeric(unlist(sio4[c(firstrow:lastrow,wagerow,surrow,taxrow), c(hhcol,surcol1,surcol2,taxcol1,taxcol2,taxcol3)])), nrow = numrows)
#final hh
fhh <- matrix(as.numeric(gsub(',','',unlist(c(sio4$H[firstrow:lastrow], "0")))), 
               nrow = numrows)
c(fhh)
#final private capital formation
fsur1 <- matrix(as.numeric(gsub(',','',unlist(c(sio4$R1[firstrow:lastrow],"0")))), 
               nrow = numrows)
c(fsur1)
fsur2 <- matrix(as.numeric(gsub(',','',unlist(c(sio4$R2[firstrow:lastrow],
                                               "0")))), 
               nrow = numrows)
c(fsur2)
fsur <- fsur1 + fsur2
View(fsur)

#final gov 
#ftax1 <- matrix(as.numeric(gsub(',','',unlist(sio4[c(firstrow:lastrow,wagerow,surrow),c(taxcol1)]))))
#ftax2 <- matrix(as.numeric(gsub(',','',unlist(sio4[c(firstrow:lastrow,wagerow,surrow),c(taxcol2)]))))
#ftax3 <- matrix(as.numeric(gsub(',','',unlist(sio4[c(firstrow:lastrow,wagerow,surrow),c(taxcol3)]))))
#ftax <- ftax1 + ftax2 + ftax3
#exports
Eval <- matrix(as.numeric(gsub(',','',unlist(c(sio4$E[firstrow:lastrow], 0)))), 
               nrow = numrows)
c(Eval)
#the only components of final demand are
fval <- fsur + Eval
c(fval)
#airexports are included

# total supply
Xval <- matrix(as.numeric(gsub(',','',unlist(c(sio4$X[firstrow:lastrow],
                                               sio4$X[wagerow])))), 
               nrow = numcols)
c(Xval)
#for ioanalysis
#fval <- rep(0,numrows)
#the matrix of intermediate transactions (see RAS)
#Zval <- matrix(as.numeric(gsub(',','',unlist(sio4[firstrow:lastcol, firstcol:lastcol]))), numint, numint)
#View(Zval)
#Zval <- matrix(c(Zval,fhh,fsur),numrows,numcols)
#Zval

```

The rows
```{r rows}

#now for the rows
W <- matrix(as.numeric(gsub(',','',
  unlist(c(sio4[wagerow, firstcol:lastcol],sio4$H[wagerow])))),
               ncol = numrows)
c(W)
S <- matrix(as.numeric(gsub(',','',
  unlist(c(sio4[surrow, firstcol:lastcol],sio4$H[surrow])))),
               ncol = numrows)
c(S)
#tax component of gva (to calculate gva when this component is exogenous)
Tv <- matrix(as.numeric(gsub(',','',
  unlist(c(sio4[taxrow, firstcol:lastcol],sio4$H[taxrow])))),
               ncol = numrows)
c(Tv)

Vval <- rbind(S, Tv)
View(Vval)

#imports
#Vval <- matrix(as.numeric(unlist(sio4[c(82), c(firstrow:numrows,hhcol)])), ncol = 71)
#for the sake of completeness, we also define
#mrow
#mhh <- matrix(as.numeric(gsub(',','',unlist(sio4[mrow,firstcol:lastcol]))))
#length(mhh)
#msur1 <- matrix(as.numeric(gsub(',','',unlist(sio4[mrow,surcol1]))))
#msur1
#msur2 <- matrix(as.numeric(gsub(',','',unlist(sio4[mrow,surcol2]))))
#msur2
#msur = msur1 + msur2
#msur
#mtax1 <- matrix(as.numeric(gsub(',','',unlist(sio4[mrow,taxcol1]))))
#mtax1
#mtax2 <- matrix(as.numeric(gsub(',','',unlist(sio4[mrow,taxcol2]))))
#mtax2
#mtax3 <- matrix(as.numeric(gsub(',','',unlist(sio4[mrow,taxcol3]))))
#mtax3
#mtax = mtax1 + mtax2 + mtax3
M <- matrix(as.numeric(gsub(',','', 
  unlist(c(sio4[mrow, firstcol:lastcol],sio4$H[mrow])))),
               ncol = numrows)
c(M)

tproductrow
Tc <- matrix(as.numeric(gsub(',','',
unlist(c(sio4[tcrow, firstcol:lastcol],sio4$H[tcrow])))),
               ncol = numrows)
c(Tc)

R <- matrix(as.numeric(gsub(',','',
  unlist(c(sio4[rrow, firstcol:lastcol],sio4$H[rrow])))),
               ncol = numrows)
c(R)
#indirect taxes (ie GST/VAT on products) and imports for final consumption
#indt <- matrix(as.numeric(sio4[c(77), c(firstrow:numrows,hhcol,77,79)]), ncol = 73)
#now for the labels (for sector labels see below)
f_label = "F" #sio4[labelrow,fcol]
f_label
E_label = "E"
E_label
M_label = "M"
M_label
V_label <- matrix(unlist(c("S", "Tv")),2,1)
V_label

sio4 <- dplyr::mutate(sio4, state = "qld")
#rm(RS_label)
RS_label <- data.table::data.table(
  sio4$state, gsub("^((\\w+\\W+){0,2}\\w+).*", "\\1", sio4$our71names))
View(RS_label)
#for HH closure
RS_label[numrows,2] <- "H"
#convert to matrix
RS_label <- matrix(unlist(RS_label[1:numrows,1:2]),
                   nrow=numrows, ncol=2)
View(RS_label)
#fV_label = sio4[c(77,81),3]
#now the inputoutput object. it seems as if RS_label has to be defined 
#within the object otherwise the dimensions clash (or perhaps this is 
#because it is the only label that needs to have the region defined 
#(hence the two items in the column))


```



RAS the transaction matrix starting from the original table.
```{r RAS 2019, include=FALSE}
sio4 <- data.table::as.data.table(readr::read_csv("sio4-19-preRAS.csv"))
View(sio4)
Zval <- matrix(as.numeric(gsub(',','',unlist(sio4[firstrow:lastrow, firstcol:lastcol]))), numint, numint)
rownames(Zval) <- sio4$our71names[firstrow:lastrow]
colnames(Zval) <- colnames(sio4)[firstcol:lastcol]
Zval[56,56]
sio3 <- data.table::as.data.table(readr::read_csv("sio3.csv"))
View(sio3)
Zvalorig <- matrix(as.numeric(gsub(',','',unlist(sio3[firstrow:lastrow, firstcol:lastcol]))), numint, numint)
Zvalorig[56,56]
for(i in c(44:71)){for 
  (j in c(44:71)) {
    Zval[i,j] = Zvalorig[i,j]
 }
}
sum(Zval[44:71,44:71])
Z = Zval
for(i in c(44:71)){for 
  (j in c(44:71)) {
    Z[i,j] = 0
 }
}
sum(Z[44:71,44:71])
Z[56,56]
Zdiff <- Zval-Z
sum(Zdiff[44:71,44:71])
Zdiff[56,56]
Udiff <- rowSums(Zdiff)
Vdiff <- colSums(Zdiff)
Udiff[56]
Vdiff[56]
U <- matrix(as.numeric(gsub(',','',unlist(sio4$U[firstrow:lastrow]))))
V <- matrix(as.numeric(gsub(',','',unlist(sio4[vrow, firstcol:lastcol]))))
#View(U-V)
U <- U-Udiff
V <- V-Vdiff
sum(Udiff)
sum(Vdiff)
U[32];V[32]
#U-V
Z <- myras(Z,U,V)
View(Z)
udiff <- U-rowSums(Z)
vdiff <- V-colSums(Z)
#colnames(udiff) <- "Ru"
#colnames(vdiff) <- "Rv"
#View(udiff)
#View(vdiff)
qldUV <- data.table::data.table(udiff, vdiff)
View(qldUV)
#write.csv(qldUV, file = './qldUV.csv',row.names = FALSE)
write.csv(Z, file = './Z.csv',row.names = FALSE)
max(abs(rowSums(Z)-U))
max(abs(colSums(Z)-V))
j <- which(abs(rowSums(Z)-U) == max(abs(rowSums(Z)-U)))
j
i <- which(abs(colSums(Z)-V) == max(abs(colSums(Z)-V)))
i

sio4$our71names[j]
sio4$our71names[i]
```



```{r Reconstuct-Z-save, include=FALSE}
# adding back in the old flows
Zfin <- Z+Zdiff
Zqld <- Zfin
View(Zqld)
#Reconstuct-Z-save, include=FALSE
#Save a time-stamped copy
Zqld19 <- Zqld
write.csv(Zqld19, file = "./Zqld19.csv", row.names = FALSE)
```
Close the model to households.
```{r typeIIa-Z)}

Zhqld <- data.frame(Zqld)
#first include the final HH demand column
Zhqld$H <- fhh[1:71]
# then the wage row
Zhqld[72,] <- wval 
Zhqld <- as.matrix(Zhqld)
View(Zhqld)

write.csv(Zhqld, file = "./Zhqld.csv", row.names = FALSE)
```


```{r for second round (if necessary)}
#U <- U - udiff
#V <- V - vdiff
sum(U + Udiff)
sum(V + Vdiff)
```




  
Here, we use the IO analysis package to create an IO object with all its essential pieces.
```{r io2019}
#install.packages("https://cran.r-project.org/package=ioanalysis")
#library(ioanalysis)
#length(fval);length(Eval);length(RS_label);length(Zhqld); length(Vval)
#the following uses Z post RAS! First remove row and column names
rownames(Zhqld) <- NULL
colnames(Zhqld) <- NULL
qld71io0 <- ioanalysis::as.inputoutput(Z=Zhqld, 
                            RS_label = RS_label,
                            f=fval, f_label=f_label,
                            E=Eval, E_label=E_label,
                            V=Vval, V_label=V_label,
                            X=Xval)
h <- ioanalysis::heatmap.io(qld71io0$Z, RS_label, sectors_x = c(1:numrows),
                sectors_y = c(1:numrows), FUN=log)
plot(h)
```

```{r}
#install.packages("ioanalysis", dependencies = TRUE)
vmult71InandOut <- ioanalysis::multipliers(qld71io19,multipliers = "output")
#wrow <- which(sio4$our71names == "W")
#vmult71InandOut <- ioanalysis::multipliers(qld71io19,multipliers = "wage", wage.row = 1)
vmult71InandOut
vmult71InandOut[24]
vmult71InandOut[28]
vmult71InandOut[44]
vmult71InandOut[46]



View(vmult71InandOut[,"output"] - colSums(qld71io19$L))

colSums(qld71io19$L)


```
Second round of RAS prep
```{r RAS prep2}
Zval <- Zhqld
Zval[56,56]
Z <- Zval
Z <- ifelse(Z==0,
               runif(1,0,1),
               Z)
#for(i in c(44:71)){for 
#  (j in c(44:71)) {
#    Z[i,j] = 0
# }
#}
Z[56,56]
Zdiff <- Zval-Z
Zdiff[56,56]
Udiff <- rowSums(Zdiff)
Vdiff <- colSums(Zdiff)
Udiff[56]
Vdiff[56]
U <- rowSums(Zh21)
V <- colSums(Zh21)
View(U-V)
U <- U-Udiff
V <- V-Vdiff
U[1];V[1]
U[72];V[72]
#check
qld71io19$A[,numrows]*xnifam21[numrows]

Zh21 <- qld71io19$A%*%diag(array(xnifam21))



U + E21 + fsur -xnifam21

Vold + Mval + R + sval + Tc + tval - t(xnifam21)

V + Mval + R + sval + Tc + tval - t(xnifam21)


sum(Mval + R + sval + Tc + tval - t(xnifam21))

sum(V)

#define new rows sums U
H21 <- qld71io19$A[,numrows]*xnifam21[numrows]
H21d <- data.table::data.table(sio4$our71names,H21, "Hdiff"= H21-fhh, "airdiff"=Eair21-Eair)
View(H21d)
H21-fhh

U21 <- xnifam21[1:numint] - H21[1:numint] - fsur[1:numint] - E21[1:numint]

U21h <- xnifam21 - H21 - fsur - E21

View(U21-rowSums(Zqld))


#(xcomp - xnifam21) - (qld71io19$A*xcomp[numrows]- H21) - (qld71io19$E - E21) - (fsur-fsur)


U <- U21h





```

```{r U check}
#check our old U makes sense
qld71io19$A[,numrows]*xcomp[numrows] - fhh

qld71io19$E - Eval

xcomp - Xval

fsur - fsur

intdiff <-  rowSums(Zhqld) - 
  (xcomp - qld71io19$A[,numrows]*xcomp[numrows] - qld71io19$E - fsur)
intdiff2 <-  rowSums(Zqld) - 
  (xcomp - Zhqld[,numrows] - qld71io19$E - fsur)

View(intdiff2)

```




```{r 2020-RAS}
#Zhqld is our base year (aka Zhqld19, but the latter is not defined)
#Zh20 is our after-shock Z, but pre-RAS
Zh20 <- qld71io19$A%*%diag(array(xnifam20))
Z <- Zh20
#Z <- ifelse(Z<1,
#               runif(1,0,1),
#               Z)
View(Z)
U <- xnifam20 - (E20 + fsur)
V <- xnifam20 - t(Mval + R + sval + Tc + tval)
View(U)
View(V)
Z <- myras(Z,U,V)
View(Z)
udiff <- U-rowSums(Z)
vdiff <- V-colSums(Z)
#colnames(udiff) <- "Ru"
#colnames(vdiff) <- "Rv"
#View(udiff)
#View(vdiff)
qldUV <- data.table::data.table(udiff, vdiff)
View(qldUV)

U <- U - udiff
V <- V - vdiff
View(U)


#check
qld71io19$A[,numrows]*xnifam21[numrows]





U + E21 + fsur -xnifam21

Vold + Mval + R + sval + Tc + tval - t(xnifam21)

V + Mval + R + sval + Tc + tval - t(xnifam21)


sum(Mval + R + sval + Tc + tval - t(xnifam21))

sum(V)

#define new rows sums U
H21 <- qld71io19$A[,numrows]*xnifam21[numrows]
H21d <- data.table::data.table(sio4$our71names,H21, "Hdiff"= H21-fhh, "airdiff"=Eair21-Eair)
View(H21d)
H21-fhh

U21 <- xnifam21[1:numint] - H21[1:numint] - fsur[1:numint] - E21[1:numint]

U21h <- xnifam21 - H21 - fsur - E21

View(U21-rowSums(Zqld))


#(xcomp - xnifam21) - (qld71io19$A*xcomp[numrows]- H21) - (qld71io19$E - E21) - (fsur-fsur)


U <- U21h





```

```{r}
Znifam20ad <- data.frame(qld71io19$A%*%diag(array(xnifam20)),
                                  E20, fsur)
Znifam20ad[73,] <- c(Mval,sio4brc$E[1], sio4brc$R1[1])
Znifam20ad[74,] <- c(R,sio4brc$E[2], sio4brc$R1[2])
Znifam20ad[75,] <- c(sval,sio4brc$E[3], sio4brc$R1[3])
Znifam20ad[76,] <- c(Tc,sio4brc$E[4], sio4brc$R1[4])
Znifam20ad[77,] <- c(tval,sio4brc$E[5], sio4brc$R1[5])
View(Znifam20ad)
Znifam20a <- as.matrix(Znifam20ad)
View(Znifam20a)


Unifam20 <- matrix(unlist(c(xnifam20[1:numrows],sio4brc$X[1:5])),ncol = 1)
View(Unifam20)

XEtot <- sum(E20) + sum(sio4brc$E) - sum(sio4brc$E[8],sio4brc$E[6])
XEtot
H20 <- qld71io19$A[,numrows]*xnifam20[numrows]
sum(H20)
XHtot <- sum(H20) + sum(sio4brc$H) - sum(sio4brc$H[8],sio4brc$H[6])
XHtot
Xsurtot <- sio4brc$R1[8] + sio4brc$R2[8]
Xsurtot
Vnifam20 <- matrix(unlist(c(xnifam20[1:numint], XHtot, XEtot, Xsurtot)),ncol = 1)
View(Vnifam20)


Z <- myras(Znifam20a,Unifam20,Vnifam20)


udiff <- Unifam20-rowSums(Z)
vdiff <- Vnifam20-colSums(Z)
qldUV <- data.table::data.table(udiff, vdiff)
View(qldUV)
View(Z)
Znifam20r <- Z
```

```{r round2}
Unifam20 <- Unifam20 - udiff
Vnifam20 <- Vnifam20 - vdiff
View(Unifam20)
View(Vnifam20)
Z <- myras(Znifam20r,Unifam20,Vnifam20)

udiff <- Unifam20-rowSums(Z)
vdiff <- Vnifam20-colSums(Z)
qldUV <- data.table::data.table(udiff, vdiff)
View(qldUV)
Znifam20rr <- Z
View(Znifam20rr)
View(Znifam20rr[,73]- Znifam20a)

```



```{r 2021 RAS prep2}
Zval <- Zhqld
Zval[56,56]
Z <- Zval
Z <- ifelse(Z==0,
               runif(1,0,1),
               Z)
#for(i in c(44:71)){for 
#  (j in c(44:71)) {
#    Z[i,j] = 0
# }
#}
Z[56,56]
Zdiff <- Zval-Z
Zdiff[56,56]
Udiff <- rowSums(Zdiff)
Vdiff <- colSums(Zdiff)
Udiff[56]
Vdiff[56]
U <- rowSums(Zh21)
V <- colSums(Zh21)
View(U-V)
U <- U-Udiff
V <- V-Vdiff
U[1];V[1]
U[72];V[72]
#check
qld71io19$A[,numrows]*xnifam21[numrows]

Zh21 <- qld71io19$A%*%diag(array(xnifam21))



U + E21 + fsur -xnifam21

Vold + Mval + R + sval + Tc + tval - t(xnifam21)

V + Mval + R + sval + Tc + tval - t(xnifam21)


sum(Mval + R + sval + Tc + tval - t(xnifam21))

sum(V)

#define new rows sums U
H21 <- qld71io19$A[,numrows]*xnifam21[numrows]
H21d <- data.table::data.table(sio4$our71names,H21, "Hdiff"= H21-fhh, "airdiff"=Eair21-Eair)
View(H21d)
H21-fhh

U21 <- xnifam21[1:numint] - H21[1:numint] - fsur[1:numint] - E21[1:numint]

U21h <- xnifam21 - H21 - fsur - E21

View(U21-rowSums(Zqld))


#(xcomp - xnifam21) - (qld71io19$A*xcomp[numrows]- H21) - (qld71io19$E - E21) - (fsur-fsur)


U <- U21h





```


Alternative RAS, where old sectors are allowed to change:
RAS prep starting from the original table.
```{r RAS 2019plus, include=FALSE}
sio4 <- data.table::as.data.table(readr::read_csv("sio4-19-preRAS.csv"))
View(sio4)
Zval <- matrix(as.numeric(gsub(',','',unlist(sio4[firstrow:lastrow, firstcol:lastcol]))), numint, numint)
rownames(Zval) <- sio4$our71names[firstrow:lastrow]
colnames(Zval) <- colnames(sio4)[firstcol:lastcol]
Zval[56,56]
sio3 <- data.table::as.data.table(readr::read_csv("sio3.csv"))
View(sio3)
Zvalorig <- matrix(as.numeric(gsub(',','',unlist(sio3[firstrow:lastrow, firstcol:lastcol]))), numint, numint)
Zvalorig[56,56]
for(i in c(44:71)){for 
  (j in c(44:71)) {
    Zval[i,j] = Zvalorig[i,j]
 }
}
sum(Zval[44:71,44:71])
Z = Zval
for(i in c(44:71)){for 
  (j in c(44:71)) {
    Z[i,j] = 0
 }
}
sum(Z[44:71,44:71])
Z[56,56]
U <- matrix(as.numeric(gsub(',','',unlist(sio4$U[firstrow:lastrow]))))
V <- matrix(as.numeric(gsub(',','',unlist(sio4[vrow, firstcol:lastcol]))))
#View(U-V)
U[32];V[32]
#U-V
Z <- myras(Z,U,V)
View(Z)
udiff <- U-rowSums(Z)
vdiff <- V-colSums(Z)
#colnames(udiff) <- "Ru"
#colnames(vdiff) <- "Rv"
#View(udiff)
#View(vdiff)
qldUV <- data.table::data.table(udiff, vdiff)
View(qldUV)
#write.csv(qldUV, file = './qldUV.csv',row.names = FALSE)
write.csv(Z, file = './Z.csv',row.names = FALSE)
max(abs(rowSums(Z)-U))
max(abs(colSums(Z)-V))
j <- which(abs(rowSums(Z)-U) == max(abs(rowSums(Z)-U)))
j
i <- which(abs(colSums(Z)-V) == max(abs(colSums(Z)-V)))
i

sio4$our71names[j]
sio4$our71names[i]
```
