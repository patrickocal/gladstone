##manual extraction of Al##
#first the benchmark status quo or complete economy
A <- glad71.IO$A
I <- diag(x=1,71)
#first complete
L <- solve(I-A)
ftot <- fbar + Ebar
xcomp <- L%*%ftot
sum(xcomp)
#check by comparing with the sum of Xbar
sum(Xbar)
#next set final demand for Al to zero
#reset A and f
A <- glad71.IO$A
ftot <- fbar + Ebar
#set final demand for Al to zero
ftot[21] <-0
#solve the model
L <- solve(I-A)
xfinAl <- L%*%ftot
#sum up the total output in this case
sum(xfinAl)
#get the difference
delxfinAl <- sum(xcomp - xfinAl)
#finally, the percentage change
perchxfinAl <- delxfinAl/sum(xcomp)*100
perchxfinAl
#difference between this and the total extraction of Al
perchxexAl - perchxfinAl


##next Al and El together##
#reset A and f
A <- glad71.IO$A
ftot <- fbar + Ebar
#set final demand for Al to zero
ftot[21] <-0
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#solve the model
L <- solve(I-A)
xfinAlEl <- L%*%ftot
#sum up the total output in this case
sum(xfinAlEl)
#get the difference
delxfinAlEl <- sum(xcomp - xfinAlEl)
#finally, the percentage change
perchxfinAlEl <- delxfinAlEl/sum(xcomp)*100
perchxfinAlEl
#difference between this and the total extraction of Al
perchxexAlEl - perchxfinAlEl

#next the extended final demand shock cases: zero for Al and, for the other sectors:
# 70 per cent extraction for El (and no exports), and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
#reset A and f
A <- glad71.IO$A
ftot <- fbar + Ebar
#set final demand for Al to zero
ftot[21] <-0
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#solve the model
L <- solve(I-A)
xfinAlEletc <- L%*%ftot
#sum up the total output in this case
sum(xfinAlEletc)
#get the difference
delxfinAlEletc <- sum(xcomp - xfinAlEletc)
#finally, the percentage change
perchxfinAlEletc <- delxfinAlEletc /sum(xcomp)*100
perchxfinAlEletc
#difference between this and the total extraction of Al
perchxexAlEletc - perchxfinAlEletc


###next we repeat extractions (all the above) assuming AlOx can replace
#its BSL contracts with exports first backward only for Al alone
#
#next set final demand for Al to zero
#just as a precaution, we reset A to the standard form
#reset A and f
A <- glad71.IO$A
ftot <- fbar + Ebar
#set final demand for Al to zero
ftot[21] <-0
#full replacement of BSL purchases of AlOx with exports
ftot[8] <-ftot[8]+219.56
#solve the model
L <- solve(I-A)
xfinAl.1 <- L%*%ftot
#sum up the total output in this case
sum(xfinAl.1)
#get the difference
delxfinAl.1 <- sum(xcomp - xfinAl.1)
#finally, the percentage change
perchxfinAl.1 <- delxfinAl.1/sum(xcomp)*100
perchxfinAl.1
#difference between this and the total extraction of Al
perchxexAl.1 - perchxfinAl.1


##next Al and El together##
#reset A and f
A <- glad71.IO$A
ftot <- fbar + Ebar
#set final demand for Al to zero
ftot[21] <-0
#full replacement of BSL purchases of AlOx with exports
ftot[8] <-ftot[8]+219.56
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#solve the model
L <- solve(I-A)
xfinAlEl.1 <- L%*%ftot
#sum up the total output in this case
sum(xfinAlEl.1)
#get the difference
delxfinAlEl.1 <- sum(xcomp - xfinAlEl.1)
#finally, the percentage change
perchxfinAlEl.1 <- delxfinAlEl.1/sum(xcomp)*100
perchxfinAlEl.1
#difference between this and the total extraction of Al
perchxexAlEl.1 - perchxfinAlEl.1

#next the extended final demand shock cases: zero for Al and, for the other sectors:
# 70 per cent extraction for El (and no exports), and 
# 10 per cent extraction for Aluminium fabrication, Transport Manufacturing, and Machinery
#(columns 22, 23 and 24 reduced to 90% of the status quo to reflect an increase in costs).
#reset A
#reset A and f
A <- glad71.IO$A
ftot <- fbar + Ebar
#set final demand for Al to zero
ftot[21] <-0
#full replacement of BSL purchases of AlOx with exports
ftot[8] <-ftot[8]+219.56
#subtract electricity exports and reduce final demand by 70%
ftot[26] <- (ftot[26] -53.4)*.3
#reduce exports of the downstream sectors by 10%
ftot[22] <- ftot[22]*.9
ftot[23] <- ftot[23]*.9
ftot[24] <- ftot[24]*.9
#solve the model
L <- solve(I-A)
xfinAlEletc.1 <- L%*%ftot
#sum up the total output in this case
sum(xfinAlEletc.1)
#get the difference
delxfinAlEletc.1 <- sum(xcomp - xfinAlEletc.1)
#finally, the percentage change
perchxfinAlEletc.1 <- delxfinAlEletc.1 /sum(xcomp)*100
perchxfinAlEletc.1
#difference between this and the total extraction of Al
perchxexAlEletc.1 - perchxfinAlEletc.1

perchxfin <- round(t(matrix(c(perchxfinAl.1,  
                       perchxfinAlEl.1,
                       perchxfinAlEletc.1,  
                       perchxfinAl,
                       perchxfinAlEl,  
                       perchxfinAlEletc),3,2)),2)
perchxfin
perchxex -perchxfin
perchxfin/perchxex
perchwex/perchxex
