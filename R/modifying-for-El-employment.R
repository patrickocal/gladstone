vecinwexAlfull <- AexAl[71,]*xexAlfull

vecinwexAlElfull <- AexAlEl[71,]*(xexAlElfull)

vecinwexAlfull[26]/Alsal
vecinwexAlElfull[26]/Alsal
  
plot(vecinwexAlfull)


vectotwexAlfull <- Acomp[71,]*xcomp - AexAl[71,]*xexAlfull

vectotwexAlElfull <- Acomp[71,]*xcomp - AexAlEl[71,]*xexAlElfull


Acomp[71.]*xcomp[26]

vectotwexAlElfull[26]

plot(vectotwexAlfull)
plot(vectotwexAlElfull)

sum(AexAlEl[71,])

#modified AlEl employee compensation
modAexAlEl71 <- 0
modAexAlEl71 <-AexAlEl[71,]*xexAlElfull/avsal
modAexAlEl71[26]
modAexAlEl71[26] <-0
modAexAlEl71 <- modAexAlEl71/sum(modAexAlEl71)*200
modAexAlEl71 <- modAexAlEl71 + AexAlEl[71,]*xexAlElfull/avsal
modAexAlEl71[26] <- AexAlEl[71,26]*xexAlElfull[26]/Alsal + 200
modAexAlEl71[26]
modAcomp71 <- Acomp[71,]*xcomp/avsal
modAcomp71[21] <- modAcomp71[21]*avsal/Alsal
modAcomp71[21]
modAcomp71[26] <- modAcomp71[26]*avsal/Alsal
modAcomp71[26]

delAexAlEl71 <- modAcomp71 -modAexAlEl71

plot(delAexAlEl71)
#modified AlElpar employee compensation
modAexAlElpar71 <- 0
modAexAlElpar71 <-AexAlElpar[71,]*xexAlElpar/avsal
modAexAlElpar71[26]
modAexAlElpar71[26] <-0
modAexAlElpar71 <- modAexAlElpar71/sum(modAexAlElpar71)*200
modAexAlElpar71 <- modAexAlElpar71 + AexAlElpar[71,]*xexAlElpar/avsal
modAexAlElpar71[26] <- AexAlElpar[71,26]*xexAlElpar[26]/Alsal + 200
modAexAlElpar71[26]
modAcomp71 <- Acomp[71,]*xcomp/avsal
modAcomp71[21] <- modAcomp71[21]*avsal/Alsal
modAcomp71[21]
modAcomp71[26] <- modAcomp71[26]*avsal/Alsal
modAcomp71[26]

delAexAlElpar71 <- modAcomp71 -modAexAlElpar71

plot(delAexAlElpar71)
#modified AlEl employee compensation
modAexAlElnone71 <- 0
modAexAlElnone71 <-AexAlElnone[71,]*xexAlElnone/avsal
modAexAlElnone71[26]
modAexAlElnone71[26] <-0
modAexAlElnone71 <- modAexAlElnone71/sum(modAexAlElnone71)*200
modAexAlElnone71 <- modAexAlElnone71 + AexAlElnone[71,]*xexAlElnone/avsal
modAexAlElnone71[26] <- AexAlElnone[71,26]*xexAlElnone[26]/Alsal + 200
modAexAlElnone71[26]
modAcomp71 <- Acomp[71,]*xcomp/avsal
modAcomp71[21] <- modAcomp71[21]*avsal/Alsal
modAcomp71[21]
modAcomp71[26] <- modAcomp71[26]*avsal/Alsal
modAcomp71[26]

delAexAlElnone71 <- modAcomp71 -modAexAlElnone71

plot(delAexAlElnone71)
