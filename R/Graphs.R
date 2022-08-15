xcomp[71]/avsal

AexAlEl[71,]*xexAlEletc
AexAlEl[71,]*xexAlEl.1
sum(xcomp)

(AexAl[71,]*xexAl.1-Acomp[71,]*xcomp)/avsal

sojlAlEl.1 <- (Acomp[71,]*xcomp- AexAlEl[71,]*xexAlEl.1)/avsal
sojlAlEl.1
sum(sojlAlEl.1)


sojlAlEl.1 <- sojlAlEl.1[-c(71)]

which.max(sojlAlEl.1)

#plot(sojlAlEl.1)
#seclab <- matrix(as.character(unlist(Glad[4:73,3])),nrow=70)
pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/job-losses.pdf", height=6, width=8)
plot(sojlAlEl.1, type='h', main='', xlab = "Sector, by number", ylab="No. Workers")
dev.off()



seclab[c(21,26,54,56,34,61,64,37)]
# Make x axis using seclab labels
axis(1, at=seq(1, 8, by=1), labels = FALSE)

par(mar = rep(2, 4))
plot(sojlAlEl[c(21,26,54,56,34,61,64,37)], xlab =seclab[c(21,26,54,56,34,61,64,37)])

