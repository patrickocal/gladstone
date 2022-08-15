rm(list=ls())
AusIO <- `_AusIO.noironmining`#read.csv("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_R/_AusIO.csv", header=FALSE)
View(AusIO)
sec <- matrix(as.numeric(unlist(AusIO[4:122,1])),nrow=119)
rows <-matrix(as.numeric(unlist(AusIO[4:122,4:125])),nrow=119)
x<-data.frame(sec,rows)
library(dplyr)
xfinr <- x %>% group_by(sec) %>% summarise_all(sum)
xragg <- xfinr[,-c(1)]
secc <- matrix(as.numeric(unlist(AusIO[1,4:125])),nrow=122)
xc <-data.frame(secc,t(xragg))
xfin <- xc %>% group_by(secc) %>% summarise_all(sum)
library(ioanalysis)
xfin <- t(xfin)
#start constructing full table
#augment columns
cregions <- matrix(rep("Glad",79),ncol=79)
xfinreg <- rbind(cregions,xfin)
#augment rows
rregions <- matrix(c("Glad",0,rep("Glad",76)))
rowlabel <- matrix(c(0,0,as.numeric(unlist(xfinr[1:76,1]))))
fulltable <-cbind(rregions,rowlabel,xfinreg)
write.csv(fulltable, file="~/OneDrive - The University of Queensland/_QTC/_Gladstone/_R/IO-70-noironmining.csv")
