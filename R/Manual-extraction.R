##manual extraction of Al##
#first the benchmark status quo or complete economy
A <- glad71.IO$A
I <- diag(x=1,71)
#first complete
L <- solve(I-A)
ftot <- fbar + Ebar +hhbar
xcomp <- L%*%ftot
sum(xcomp)
#next backward only
#reset A
A <- glad71.IO$A
A[,21] <- rep(0,71)
L <- solve(I-A)
ftot <- fbar + Ebar +hhbar
ftot[21] <-0
xexinAl <- L%*%ftot
sum(xexinAl)
delexinAl <- sum(xcomp - xexinAl)
delexinAl/sum(xcomp)
#next total
#reset A and f
A <- glad71.IO$A
ftot <-fbar + Ebar + hhbar
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
L <- solve(I-A)
ftot <- fbar + Ebar +hhbar
ftot[21] <-0
xextotAl <- L%*%ftot
sum(xextotAl)
delextotAl <- sum(xcomp - xextotAl)
delextotAl
delextotAl/sum(xcomp)
(delextotAl - delexinAl)/sum(xcomp)

##next Al and El together##
#reset A and f
A <- glad71.IO$A
ftot <-fbar + Ebar + hhbar
#backward only
A[,21] <- rep(0,71)
A[,26] <- rep(0,71)
#A[21,] <- rep(0,71)
L <- solve(I-A)
ftot <- fbar + Ebar +hhbar
ftot[21] <-0
ftot[26] <-0
xexinAlEl <- L%*%ftot
sum(xexinAlEl)
delexinAlEl <- sum(xcomp - xexinAlEl)
delexinAlEl/sum(xcomp)
#next total for Al and inputs/backward linkage only for El
#reset A and f
A <- glad71.IO$A
ftot <-fbar + Ebar + hhbar
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#only inputs (backward linkages) for electricity
A[,26] <- rep(0,71)
L <- solve(I-A)
ftot <- fbar + Ebar +hhbar
ftot[21] <-0
ftot[26] <-0
xextotAlinEl <- L%*%ftot
sum(xextotAlinEl)
delextotAlinEl <- sum(xcomp - xextotAlinEl)
delextotAlinEl/sum(xcomp)
(delextotAlinEl - delexinAlEl)/sum(xcomp)

###next we repeat extractions assuming AlOx can replace its BSL contracts with exports
#first backward only for Al alone
#reset A
A <- glad71.IO$A
A[,21] <- rep(0,71)
L <- solve(I-A)
ftot <- fbar + Ebar +hhbar
ftot[21] <-0
ftot[8] <-ftot[8]+219.56
xexinAl.1 <- L%*%ftot
sum(xexinAl.1)
delexinAl.1 <- sum(xcomp - xexinAl.1)
delexinAl.1/sum(xcomp)
#next total
#reset A and f
A <- glad71.IO$A
ftot <-fbar + Ebar + hhbar
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
L <- solve(I-A)
ftot <- fbar + Ebar +hhbar
ftot[21] <-0
ftot[8] <-ftot[8]+219.56
xextotAl.1 <- L%*%ftot
sum(xextotAl.1)
delextotAl.1 <- sum(xcomp - xextotAl.1)
delextotAl.1/sum(xcomp)
(delextotAl.1 - delexinAl.1)/sum(xcomp)

##next Al and El together##
#reset A and f
A <- glad71.IO$A
ftot <-fbar + Ebar + hhbar
#backward only
A[,21] <- rep(0,71)
A[,26] <- rep(0,71)
#A[21,] <- rep(0,71)
L <- solve(I-A)
ftot <- fbar + Ebar +hhbar
ftot[21] <-0
ftot[26] <-0
ftot[8] <-ftot[8]+219.56
xexinAlEl.1 <- L%*%ftot
sum(xexinAlEl.1)
delexinAlEl.1 <- sum(xcomp - xexinAlEl.1)
delexinAlEl.1/sum(xcomp)
#next total for Al and inputs/backward linkage only for El
#reset A and f
A <- glad71.IO$A
ftot <-fbar + Ebar + hhbar
#total extraction of Al
A[,21] <- rep(0,71)
A[21,] <- rep(0,71)
#only inputs (backward linkages) for electricity
A[,26] <- rep(0,71)
L <- solve(I-A)
ftot <- fbar + Ebar +hhbar
ftot[21] <-0
ftot[26] <-0
ftot[8] <-ftot[8]+219.56
xextotAlinEl.1 <- L%*%ftot
sum(xextotAlinEl.1)
delextotAlinEl.1 <- sum(xcomp - xextotAlinEl.1)
delextotAlinEl.1/sum(xcomp)
(delextotAlinEl.1 - delexinAlEl.1)/sum(xcomp)
#effect of AlOx export response
delexinAl - delexinAl.1
delextotAl - delextotAl.1
delexinAlEl - delexinAlEl.1
delextotAlinEl - delextotAlinEl.1

(delexinAl - delexinAl.1)/delexinAl
(delextotAl - delextotAl.1)/delextotAl
(delexinAlEl - delexinAlEl.1)/delexinAlEl
(delextotAlinEl - delextotAlinEl.1)/delextotAlinEl

(delexinAl - delexinAl.1)/sum(xcomp)
(delextotAl - delextotAl.1)/sum(xcomp)
(delexinAlEl - delexinAlEl.1)/sum(xcomp)
(delextotAlinEl - delextotAlinEl.1)/sum(xcomp)

maxdel <- max(delexinAl, delexinAl, delexinAl.1, delextotAl, 
              delextotAl.1, delexinAlEl, delexinAlEl.1,delextotAlinEl, delextotAlinEl.1)
mindel <- min(delexinAl, delexinAl, delexinAl.1, delextotAl, 
              delextotAl.1, delexinAlEl, delexinAlEl.1,delextotAlinEl, delextotAlinEl.1)
maxdel/sum(xcomp)
mindel/sum(xcomp)

(xcomp[71] - xexinAl.1[71])/xcomp[71]

(xcomp[71] - xextotAlinEl.1[71])/xcomp[71]
xextotAl.1[71]-sum(wbar)
xextotAlinEl[71]-sum(wbar)
t(xextotAlinEl)
sum(wbar)
