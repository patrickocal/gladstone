##extraction of sectors##
#first the benchmark status quo or complete economy
#let us define the matrix of technical coefficients
A <- vglad71.IO$A
#and final output (recall that final HH is part of the intermediate sector) 
ftot <- Eval
I <- diag(x=1,73)
L <- solve(I-A)
xcomp <- L%*%ftot
sum(xcomp)
#check by comparing with the sum of Xval
sum(Xval)
xcomp[71]
xcomp -A%*%xcomp-Eval
#therefore
wcomp <- A[71,]*xcomp
xcomp[71]
sum(wcomp)
wshare <- sum(wcomp)/sum(xcomp)
wshare
scomp <- A[72,]*xcomp
xcomp[72]
sum(scomp)
sshare <- sum(scomp)/sum(xcomp)
sshare
tcomp <- A[73,]*xcomp
sum(tcomp)
tshare <- sum(tcomp)/sum(xcomp)
tshare
gvacomp <- wcomp + scomp + tcomp 
sum(gvacomp)
sum(xcomp[71:73])
gvashare <- sum(gvacomp)/sum(xcomp)
gvashare
wshare/gvashare
#next total extraction of Al
#reset A
A <- vglad71.IO$A
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#let's view the relevant part of our new intermediate matrix A using a heatmap
obj = matrix(c(A,ftot),ncol=dim(A)[1]+1)
hexAl <- heatmap.io(obj, RS_label, sectors_x = c(5:34),
                    sectors_y = c(5:34), FUN=log, max = 7)
plot(hexAl)
#reset final demand
ftot <- Eval
#now let's set final demand for Al to zero (so everything comes from imports )
ftot[21] <-0
#now let's solve the model
L <- solve(I-A)
xexAl <- L%*%ftot
#therefore
wexAl <- A[71,]*xexAl
xexAl[71]
sum(wexAl)
wshare <- sum(wexAl)/sum(xexAl)
wshare
sexAl <- A[72,]*xexAl
xexAl[72]
sum(sexAl)
sshareexAl <- sum(sexAl)/sum(xexAl)
sshareexAl
texAl <- A[73,]*xexAl
sum(texAl)
tshareexAl <- sum(texAl)/sum(xexAl)
tshareexAl
gvaexAl <- wexAl + sexAl + texAl 
sum(gvaexAl)
sum(xexAl[71:73])
gvashareexAl <- sum(gvaexAl)/sum(xexAl)
gvashareexAl
wshareexAl/gvashareexAl
#sum up the new output vector
sum(xexAl)
#compute the difference
delexAl <- sum(xcomp - xexAl)
#compute and then display the percentage loss
perchxexAl <- round(delexAl/sum(xcomp)*100,2)
perchxexAl
#this completes the extraction of Al
#little check
xexAl -A%*%xexAl-ftot
#indirect impact on value added
A[71,]*(xcomp-xexAl)
sum(A[71,]*(xcomp-xexAl))
#the following includes direct and indirect
xcomp[71]-xexAl[71]
(xcomp[71]-xexAl[71])/xcomp[71]






#
##next Al and El together##
#reset A
A <- vglad71.IO$A
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#total extraction of El
A[,26] <- A[,26]*rep(0.3,71)
A[26,] <- A[26,]*rep(0.3,71)
#reset ftot
ftot <- Eval
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
#little check
xexAlEl -A%*%xexAlEl-ftot
#comparison with just Al
perchxexAlEl - perchxexAl

#next the extended extraction cases: total for Al and, for the other sectors:
# 70 per cent extraction for El, and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
A <- vglad71.IO$A
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
ftot <- Eval
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

####next we repeat all three extractions assuming AlOx can replace
#its BSL contracts with exports####
#total extraction of Al
#reset A
A <- vglad71.IO$A
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#let's view the relevant part of our new intermediate matrix A using a heatmap
obj = matrix(c(A,ftot),ncol=dim(A)[1]+1)
hexAl <- heatmap.io(obj, RS_label, sectors_x = c(5:34),
                    sectors_y = c(7:34), FUN=log, max = 7)
plot(hexAl)
#reset final demand
ftot <- Eval
#now let's set final demand for Al to zero (so everything comes from imports )
ftot[21] <-0
#full replacement of BSL purchases of AlOx with exports
ftot[8] <-ftot[8]+219.56
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

##next Al and El together##
#reset A
A <- vglad71.IO$A
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#total extraction of El
A[,26] <- A[,26]*rep(0.3,71)
A[26,] <- A[26,]*rep(0.3,71)
#reset ftot
ftot <- Eval
ftot[21] <-0
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#full replacement of BSL purchases of AlOx with exports
ftot[8] <-ftot[8]+219.56
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
A <- vglad71.IO$A
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
ftot <- Eval
ftot[21] <-0
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#full replacement of BSL purchases of AlOx with exports
ftot[8] <-ftot[8]+219.56
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

vperchxex <- t(matrix(c(perchxexAl.1,  
                       perchxexAlEl.1,
                       perchxexAlEletc.1,  
                       perchxexAl,
                       perchxexAlEl,  
                       perchxexAlEletc),3,2))

vperchxex



perchvexAl.1 <- (xcomp[71] - xexAl.1[71])/xcomp[71]*100
perchvexAlEl.1 <- (xcomp[71] - xexAlEl.1[71])/xcomp[71]*100
perchvexAlEletc.1 <- (xcomp[71] - xexAlEletc.1[71])/xcomp[71]*100

perchvexAl <- (xcomp[71]-xexAl[71])/xcomp[71]*100
perchvexAlEl <- (xcomp[71]-xexAlEl[71])/xcomp[71]*100
perchvexAlEletc <-(xcomp[71] - xexAlEletc[71])/xcomp[71]*100

perchvex <- t(matrix(c(perchvexAl.1,  
                       perchvexAlEl.1,
                       perchvexAlEletc.1,  
                       perchvexAl,
                       perchvexAlEl,  
                       perchvexAlEletc),3,2))
perchvex
perchvexAl.1/perchxexAl.1

round(perchvex/vperchxex,1)

perchwex/vperchxex
perchwex/perchxex
perchwex/perchvex

dperchwex <- matrix (rep(0,12), 3,4)
dperchwex



