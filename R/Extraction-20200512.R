##extraction of sectors##
#first the benchmark status quo or complete economy
#let us define the matrix of technical coefficients
A <- glad71.IO$A
#and final output (recall that final HH is part of the intermediate sector) 
ftot <- fbar + Ebar
I <- diag(x=1,71)
L <- solve(I-A)
xcomp <- L%*%ftot
sum(xcomp)
#check by comparing with the sum of Xbar
sum(Xbar)
xcomp[71]
xcomp -A%*%xcomp-ftot
#therefore
A[71,]*xcomp
sum(A[71,]*xcomp)
###total extraction of Al###
#reset A
A <- glad71.IO$A
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#let's view the relevant part of our new intermediate matrix A using a heatmap
obj = matrix(c(A,ftot),ncol=dim(A)[1]+1)
hexAl <- heatmap.io(obj, RS_label, sectors_x = c(5:34),
                     sectors_y = c(5:34), FUN=log, max = 7)
plot(hexAl)
#reset final demand
ftot <- fbar + Ebar
#now let's set final demand for Al to zero (so everything comes from imports )
ftot[21] <-0
#now let's solve the model
L <- solve(I-A)
xexAl <- L%*%ftot
#sum up the new output vector
sum(xexAl)
#compute the difference
delexAl <- sum(xcomp - xexAl)
#compute and then display the percentage loss
perchxexAl <- round(delexAl/sum(xcomp)*100,2)
perchxexAl
#this completes the extraction of Al
#indirect impact on compensation of employees
A[71,]*(xcomp-xexAl)
sum(A[71,]*(xcomp-xexAl))
xcomp[71]-xexAl[71]

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
#reset f
ftot <- fbar + Ebar
ftot[21] <-0
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#solve the model
L <- solve(I-A)
xexAletc <- L%*%ftot
#sum up the new output vector
sum(xexAletc)
#compute the difference
delexAletc <- sum(xcomp - xexAletc)
#compute and then display the percentage loss
perchxexAletc <- round(delexAletc/sum(xcomp)*100,2)
perchxexAletc
#this completes the extraction of Al, El and the other sectors
#comparison with Al ie increment over that extraction
perchxexAletc - perchxexAl



##next Al and El together##
#reset A
A <- glad71.IO$A
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#total extraction of El
A[,26] <- A[,26]*rep(0.3,71)
A[26,] <- A[26,]*rep(0.3,71)
#reset ftot
ftot <- fbar + Ebar
ftot[21] <-0
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#solve the model
L <- solve(I-A)
xexAlEl <- L%*%ftot
#sum up the new output vector
sum(xexAlEl)
#compute the difference
delexAlEl <- sum(xcomp - xexAlEl)
#compute and then display the percentage loss
perchxexAlEl <- round(delexAlEl/sum(xcomp)*100,2)
perchxexAlEl
#this completes the extraction of Al and El
#comparison with just Al
perchxexAlEl - perchxexAl

#next the extended extraction cases: total for Al and, for the other sectors:
# 70 per cent extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- glad71.IO$A
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
# 70 percent total extraction of El
A[,26] <- A[,26]*.3
A[26,] <- A[26,]*.3
# 10 percent total extraction of the three potential downstream sectors
A[,22] <- .9*A[,22]
A[,23] <- .9*A[,23]
A[,24] <- .9*A[,24]
#reset f
ftot <- fbar + Ebar
ftot[21] <-0
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#solve the model
L <- solve(I-A)
xexAlEletc <- L%*%ftot
#sum up the new output vector
sum(xexAlEletc)
#compute the difference
delexAlEletc <- sum(xcomp - xexAlEletc)
#compute and then display the percentage loss
perchxexAlEletc <- round(delexAlEletc/sum(xcomp)*100,2)
perchxexAlEletc
#this completes the extraction of Al, El and the other sectors
#comparison with Al + El ie increment over that extraction
perchxexAlEletc - perchxexAlEl

####next we repeat all three extractions assuming AlOx and El can replace
#its BSL contracts with exports####
#total extraction of Al
#reset A
A <- glad71.IO$A
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#let's view the relevant part of our new intermediate matrix A using a heatmap
obj = matrix(c(A,ftot),ncol=dim(A)[1]+1)
hexAl <- heatmap.io(obj, RS_label, sectors_x = c(5:34),
                     sectors_y = c(7:34), FUN=log, max = 7)
plot(hexAl)
#reset final demand
ftot <- fbar + Ebar
#now let's set final demand for Al to zero (so everything comes from imports )
ftot[21] <-0
#full replacement of BSL purchases of AlOx with exports
Zbar[8,21]
ftot[8] <-ftot[8]+Zbar[8,21]
#and full replacement of BSL purchases of electricity with exports
Zbar[26,21]
ftot[26] <- ftot[26] + Zbar[26,21]
#now let's solve the model
L <- solve(I-A)
xexAl.1 <- L%*%ftot
#sum up the new output vector
sum(xexAl.1)
#compute the difference
delexAl.1 <- sum(xcomp - xexAl.1)
#compute and then display the percentage loss
perchxexAl.1 <- round(delexAl.1/sum(xcomp)*100,2)
perchxexAl.1
#this completes the extraction of Al
#comparison with Al alone
perchxexAl - perchxexAl.1
#indirect impact on compensation of employees
A[71,]*(xcomp-xexAl.1)
sum(A[71,]*(xcomp-xexAl.1))
xcomp[71]-xexAl.1[71]


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
#reset f
ftot <- fbar + Ebar
ftot[21] <-0
#full replacement of BSL purchases of AlOx with exports
Zbar[8,21]
ftot[8] <-ftot[8]+Zbar[8,21]
#and full replacement of BSL purchases of electricity with exports
Zbar[26,21]
ftot[26] <- ftot[26] + Zbar[26,21]
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#solve the model
L <- solve(I-A)
xexAletc.1 <- L%*%ftot
#sum up the new output vector
sum(xexAletc.1)
#compute the difference
delexAletc.1 <- sum(xcomp - xexAletc.1)
#compute and then display the percentage loss
perchxexAletc.1 <- round(delexAletc.1/sum(xcomp)*100,2)
perchxexAletc.1
#this completes the extraction of Al, El and the other sectors
#comparison with Al ie increment over that extraction
perchxexAletc.1 - perchxexAl.1


###next Al and El together###
#reset A
A <- glad71.IO$A
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#total extraction of El
A[,26] <- A[,26]*rep(0.3,71)
A[26,] <- A[26,]*rep(0.3,71)
#reset ftot
ftot <- fbar + Ebar
ftot[21] <-0
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#full replacement of BSL purchases of AlOx with exports
Zbar[8,21]
ftot[8] <-ftot[8]+Zbar[8,21]
#solve the model
L <- solve(I-A)
xexAlEl.1 <- L%*%ftot
#sum up the new output vector
sum(xexAlEl.1)
#compute the difference
delexAlEl.1 <- sum(xcomp - xexAlEl.1)
#compute and then display the percentage loss
perchxexAlEl.1 <- round(delexAlEl.1/sum(xcomp)*100,2)
perchxexAlEl.1
#this completes the extraction of Al and El
#comparison with just Al
perchxexAlEl - perchxexAlEl.1
#comparison with Al.1
perchxexAlEl.1 - perchxexAl.1
perchxexAlEl - perchxexAl
#next the extended extraction cases: total for Al and, for the other sectors:
# 70 per cent extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- glad71.IO$A
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
# 70 percent total extraction of El
A[,26] <- A[,26]*.3
A[26,] <- A[26,]*.3
# 10 percent total extraction of the three potential downstream sectors
A[,22] <- .9*A[,22]
A[,23] <- .9*A[,23]
A[,24] <- .9*A[,24]
#reset f
ftot <- fbar + Ebar
ftot[21] <-0
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#full replacement of BSL purchases of AlOx with exports
Zbar[8,21]
ftot[8] <-ftot[8]+Zbar[8,21]
#solve the model
L <- solve(I-A)
xexAlEletc.1 <- L%*%ftot
#sum up the new output vector
sum(xexAlEletc.1)
#compute the difference
delexAlEletc.1 <- sum(xcomp - xexAlEletc.1)
#compute and then display the percentage loss
perchxexAlEletc.1 <- round(delexAlEletc.1/sum(xcomp)*100,2)
perchxexAlEletc.1
#this completes the extraction of Al, El and the other sectors
#comparison 
perchxexAlEletc - perchxexAlEletc.1
perchxexAlEletc.1 - perchxexAl.1
perchxexAlEletc - perchxexAl

#full comparison
perchxex <- round(t(matrix(c(perchxexAl.1,
                       perchxexAletc.1,
                       perchxexAlEl.1,
                       perchxexAlEletc.1,  
                       perchxexAl,
                       perchxexAletc,
                       perchxexAlEl,  
                       perchxexAlEletc),4,2)),1)
perchxex
perchxeximp <- round(t(matrix(c(perchxexAl.1,
                             perchxexAlEl.1,
                             perchxexAl,
                             perchxexAlEl),2,2)),1)
perchxeximp
perchxexparimp <- round(t(matrix(c(perchxexAletc.1,
perchxexAlEletc.1,  
perchxexAletc,
perchxexAlEletc),2,2)),1)
perchxexparimp

(xcomp[71] - xexAl.1[71])/xcomp[71]

(xcomp[71] - xexAlEl.1[71])/xcomp[71]
(xcomp[71]-xexAlEl[71])/xcomp[71]
(xcomp[71]-xexAl[71])/xcomp[71]



perchwexAl.1 <- (xcomp[71] - xexAl.1[71])/xcomp[71]*100
perchwexAletc.1 <-(xcomp[71] - xexAletc.1[71])/xcomp[71]*100
perchwexAlEl.1 <- (xcomp[71] - xexAlEl.1[71])/xcomp[71]*100
perchwexAlEletc.1 <- (xcomp[71] - xexAlEletc.1[71])/xcomp[71]*100

perchwexAl <- (xcomp[71]-xexAl[71])/xcomp[71]*100
perchwexAletc <-(xcomp[71] - xexAletc[71])/xcomp[71]*100
perchwexAlEl <- (xcomp[71]-xexAlEl[71])/xcomp[71]*100
perchwexAlEletc <-(xcomp[71] - xexAlEletc[71])/xcomp[71]*100


perchwex <- round(t(matrix(c(perchwexAl.1, 
                             perchwexAletc.1,
                             perchwexAlEl.1,
                             perchwexAlEletc.1,
                             perchwexAl,
                             perchwexAletc,
                             perchwexAlEl,
                             perchwexAlEletc),4,2)),1)
perchwex
perchwexAl.1/perchxexAl.1

perchwexAlEl/perchxexAlEl
round(perchwex/perchxex*100,1)
totjobex <- perchwex/100*30000
totjobex
dirjobex <- matrix(c(960, 960, 1060, 1060, 1260, 1260, 1360, 1360),2,4)
dirjobex
indirjobex <- totjobex -dirjobex
indirjobex

indirjobexAl <- perch/100*30000 -910
totaljobexAl.1 <- indirjoblossexAl.1 + dirjoblossexAl.1
joblossexAl.1 <- matrix(c(dirjoblossexAl.1,indirjoblossexAl.1, totaljoblossexAl.1),1,3)
jobloss

perchwex
#so 378 indirect job losses in the canonical best case scenario, where only Al is extracted
#and exports replace direct transaction losses

rlabel <- matrix(as.character(c("Upstream exports respond","Upstream exports constant")))
clabel <- matrix(as.character(c("Al", "Al", "Al,El",  "Al,El,etc"))) 


dperchwex <- as.data.frame(perchwex,row.names = rlabel, col.names =clabel)

rownames(perchwex) <- rlabel
colnames(perchwex) <- clabel

perchwex

xtable(perchwex,caption = "Percentage change in Total Employee compensation for each scenario",
       align = "||l|c|c|c|c||")  




c2label <- matrix(as.character(c("Only Al", "Both Al and El")))
rownames(perchxeximp) <- rlabel
colnames(perchxeximp) <- c2label
res.table <-xtable(perchxeximp,caption = "Percentage change in Total Output for each scenario",
       align = "lrr",digits=1)
print(res.table, scalebox=0.9, caption.placement = "top")




