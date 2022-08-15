seclab <- Glad[4:73,3]
seclab <- matrix(as.character(unlist(Glad[4:73,3])),ncol=70)
seclab

sojlAlfull <- (Acomp[71,]*xcomp- AexAl[71,]*xexAlfull)/avsal
sojlAlfull[21]
sojlAlfull[21] <- 920
sojlAlfull[26] <- sojlAlfull[26]*avsal/Alsal
sojlAlfull[26]
length(sojlAlfull)
sojlAlfull <- sojlAlfull[-c(73)]
sojlAlfull <- sojlAlfull[-c(72)]
sojlAlfull <- sojlAlfull[-c(71)]
length(sojlAlfull)
sojlAlfull

sojlAlfull <- (Acomp[71,]*xcomp- AexAl[71,]*xexAlfull)/avsal
sojlAlfull <- sojlAlfull[-c(73)]
sojlAlfull <- sojlAlfull[-c(72)]
sojlAlfull <- sojlAlfull[-c(71)]

#modified AlEl employee compensation
modAexAlEl71 <- rep(0,73)
modAexAlEl71 <-AexAlEl[71,]*xexAlElfull/avsal + modAexAlEl71
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
sojlAlElfull <- modAcomp71 -modAexAlEl71
length(sojlAlElfull)
sojlAlElfull <- sojlAlElfull[-c(73)]
sojlAlElfull <- sojlAlElfull[-c(72)]
sojlAlElfull <- sojlAlElfull[-c(71)]
length(sojlAlfull)
sojlAlElfull



rm(sojlAlpar)
sojlAlpar <- (Acomp[71,]*xcomp- AexAlpar[71,]*xexAlpar)/avsal
sojlAlpar[21] <- 920
sojlAlpar[26] <- sojlAlpar[26]*avsal/Alsal
sojlAlpar[26]
length(sojlAlpar)
sojlAlpar <- sojlAlpar[-c(73)]
sojlAlpar <- sojlAlpar[-c(72)]
sojlAlpar <- sojlAlpar[-c(71)]
length(sojlAlpar)
sojlAlpar

#modified AlEl employee compensation
rm(modAexAl71)
modAexAl71 <-AexAlnone[71,]*xexAlnone/avsal
modAexAl71[26]
modAexAl71[26] <-0
modAexAl71 <- modAexAl71/sum(modAexAl71)*91
modAexAl71 <- modAexAl71 + AexAl[71,]*xexAlnone/avsal
modAexAl71[26] <- AexAl[71,26]*xexAlnone[26]/Alsal + 91
modAexAl71[26]
modAcomp71 <- Acomp[71,]*xcomp/avsal
modAcomp71[21] <- modAcomp71[21]*avsal/Alsal
modAcomp71[21]
modAcomp71[26] <- modAcomp71[26]*avsal/Alsal
modAcomp71[26]
sojlAlnone <- modAcomp71 -modAexAl71
sojlAlnone[26]
length(sojlAlnone)
sojlAlnone <- sojlAlnone[-c(73)]
sojlAlnone <- sojlAlnone[-c(72)]
sojlAlnone <- sojlAlnone[-c(71)]
length(sojlAlnone)
sojlAlnone

rm(sojlAlnone)
sojlAlnone <- (Acomp[71,]*xcomp- AexAlnone[71,]*xexAlnone)/avsal
sojlAlnone[21] <- 920
sojlAlnone[26] <- sojlAlnone[26]*avsal/Alsal
#since this too is too high
sojlAlnone[26] <- sojlAlElfull[26]
sojlAlnone[26]
length(sojlAlnone)
sojlAlnone <- sojlAlnone[-c(73)]
sojlAlnone <- sojlAlnone[-c(72)]
sojlAlnone <- sojlAlnone[-c(71)]
length(sojlAlnone)
sojlAlnone


#modified AlElpar employee compensation
rm(modAexAlElpar71)
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
sojlAlElpar <- modAcomp71 -modAexAlElpar71
sojlAlElpar[26]
length(sojlAlElpar)
sojlAlElpar <- sojlAlElpar[-c(73)]
sojlAlElpar <- sojlAlElpar[-c(72)]
sojlAlElpar <- sojlAlElpar[-c(71)]
length(sojlAlElpar)
sojlAlElpar

#modified AlEl employee compensation
rm(modAexAlElnone71)
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
sojlAlElnone <- modAcomp71 -modAexAlElnone71
length(sojlAlElnone)
sojlAlElnone <- sojlAlElnone[-c(73)]
sojlAlElnone <- sojlAlElnone[-c(72)]
sojlAlElnone <- sojlAlElnone[-c(71)]
length(sojlAlElnone)
sojlAlElnone

library(ggpubr)

require(gridExtra)


#some minor corrections
#sojlAlfull[58] <- 0
#sojlAlnone[8] <- sojlAlElnone[8]
#sojlAlElpar[8] <- 0
#sojlAlElfull[8] <- 0

df5 <- data.frame(Scenario=rep(c("1a", "2a","1b","2b","1c","2c"), each=70),
                  Sectors=rep(1:70,6),
                  Joblosses=c(sojlAlfull,sojlAlElfull,sojlAlpar,sojlAlElpar,sojlAlnone,sojlAlElnone))
df1 <- data.frame(Sectors=rep(1:70,1),
                  Joblosses=c(sojlAlfull))
print(df1)
df4 <- data.frame(Sectors=rep(1:70,1),
                  Joblosses=c(sojlAlElnone))
print(df4)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scale_colour_manual(values=cbbPalette)

p = ggplot(df5,aes(Sectors,Joblosses),xlim=.5:70.5) +
  geom_segment(aes(x=Sectors,xend=Sectors,y=0,yend=Joblosses),color="black",size=1/10,alpha=1/10) +
  geom_point(aes(shape=Scenario),size=1.3)+
  geom_point(aes(color=Scenario),size=.8)+
  geom_line(data = df4, color = "#EE82EE",size=7/10)+
  geom_line(data = df1, color = "#FF8000",size=7/10)+
  scale_colour_hue(l=60)+
  scale_fill_manual(values=cbbPalette)+
  scale_x_continuous(limits = c(.85,70.15), expand = c(0, 0)) +
  guides(color = guide_legend(nrow=1),
         shape = guide_legend(nrow=1))
#  theme_bw()
#  theme(panel.grid = element_blank(),
#        panel.border = element_blank())
plot(p)
ggpar(p,
      #x.text.angle = 60
      font.family="Times",legend="top",font.legend=14,legend.box="vertical")
