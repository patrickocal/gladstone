sojlAl.1 <- (Acomp[71,]*xcomp- AexAl[71,]*xexAl.1)/avsal
sojlAl.1
sum(sojlAl.1)
sojlAl.1 <- sojlAl.1[-c(71)]
dfAl.1 <- data.frame(x = seq(1, 70, by = 1),
                  y = sojlAl.1)



sojlAlpar <- (Acomp[71,]*xcomp- AexAlpar[71,]*xexAl.par.etc)/avsal
sojlAlpar
sum(sojlAlpar)
sojlAlpar <- sojlAlpar[-c(71)]
dfAlpar <- data.frame(x = seq(1, 70, by = 1),
                  y = sojlAlpar)

sojlAlnone <- (Acomp[71,]*xcomp- AexAlnone[71,]*xexAletc)/avsal
sojlAlnone
sum(sojlAlnone)
sojlAlnone <- sojlAlnone[-c(71)]
dfAlnone <- data.frame(x = seq(1, 70, by = 1),
                      y = sojlAlnone)

sojlAlEl.1 <- (Acomp[71,]*xcomp- AexAlEl[71,]*xexAlEl.1)/avsal
sojlAlEl.1
sum(sojlAlEl.1)
sojlAlEl.1 <- sojlAlEl.1[-c(71)]
dfAlEl.1 <- data.frame(x = seq(1, 70, by = 1),
                     y = sojlAlEl.1)


sojlAlElpar <- (Acomp[71,]*xcomp- AexAlElpar[71,]*xexAlEl.par.etc)/avsal
sojlAlElpar
sum(sojlAlElpar)
sojlAlElpar <- sojlAlElpar[-c(71)]
dfAlElpar <- data.frame(x = seq(1, 70, by = 1),
                      y = sojlAlElpar)

sojlAlElnone <- (Acomp[71,]*xcomp- AexAlElnone[71,]*xexAlEletc)/avsal
sojlAlElnone
sum(sojlAlElnone)
sojlAlElnone <- sojlAlElnone[-c(71)]
dfAlElnone <- data.frame(x = seq(1, 70, by = 1),
                      y = sojlAlElnone)


library(ggpubr)

#now for the barplot
jlAl.1 <- matrix(sojlAl.1,ncol=70)
seclab<- matrix(as.character(seclab))
colnames(jlAl.1) <- seclab
jlAl.1 <- round(jlAl.1,0)


#matrix of bar heights for the bar plot
jlall <- matrix(c(sojlAl.1,sojlAlEl.1,sojlAlnone,sojlAlElnone),ncol=70)

end_point = 0.5 + ncol(jlall) + ncol(jlall) - 1 #this is the line which does the trick (together with barplot "space = 1" parameter)
colnames(jlall) <- seclab

barplot(jlAl.1,
        col = "cyan", 
        main = "",
        ylab = "Total Job Losses", ylim = c(0,5 + max(jlall)),
        names.arg=NULL,
        space = 1)
#rotate 60 degrees (srt = 60)
text(seq(1.5, end_point, by = 2), par("usr")[3]-0.25, 
     srt = 60, adj = 1, xpd = TRUE,
     labels = paste(seclab), cex = 0.35)

plot(jlall,beside=TRUE,type='l')


barplot(height = jlall,
        space = 1,
        names.arg=seclab , 
        col='cyan', 
        )
install.packages('ggbarplot')
library(ggbarplot)
# Plot with multiple groups
# +++++++++++++++++++++

# Create some data
df5 <- data.frame(Scenario=rep(c("1a", "2a","1c","2c"), each=70),
                  sec=rep(1:70,4),
                  Joblosses=c(sojlAl.1,sojlAlEl.1,sojlAlnone,sojlAlElnone))
print(df5)
#>   supp dose  len
#> 1   VC D0.5  6.8
#> 2   VC   D1 15.0
#> 3   VC   D2 33.0
#> 4   OJ D0.5  4.2
#> 5   OJ   D1 10.0
#> 6   OJ   D2 29.5

# Change position: Interleaved (dodged) bar plot
jlgg <- ggbarplot(df4, "sec", "Joblosses",
          fill = "Scenario", color = "Scenario", palette = c("blue","red","green","magenta"),
          position = position_dodge(.7),
          size=.1
          )
ggpar(jlgg,x.text.angle = 60,font.x=1,font.family="Times")
plot(jlgg)
pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/jlgg.pdf", height=6, width=8)
dev.off()
ggdotchart(df2, x = "sec", y = "Joblosses",
           color = "Scenario",                                # Color by groups
           palette = c("#00AFBB", "#E7B800", "#FC4E07","blue"), # Custom color palette
           sorting = "ascending",                        # Sort value in descending order
           add = "segments", data = seclab,
           position = position_dodge(.7),
           ggtheme = theme_pubr()                        # ggplot2 theme
)


seclab




text(seq(1.5, end_point, by = 2), par("usr")[3]-0.25, 
     srt = 60, adj = 1, xpd = TRUE,
     labels = paste(seclab), cex = 0.35)





end_point = 0.5 + ncol(jlall) + ncol(jlall) - 1 #this is the line which does the trick (together with barplot "space = 1" parameter)
colnames(jlall) <- seclab

barplot(jlAl.1,
        col = "cyan", 
        main = "",
        ylab = "Total Job Losses", ylim = c(0,5 + max(jlall)),
        names.arg=NULL,
        space = 1)
#rotate 60 degrees (srt = 60)
text(seq(1.5, end_point, by = 2), par("usr")[3]-0.25, 
     srt = 60, adj = 1, xpd = TRUE,
     labels = paste(seclab), cex = 0.35)






ggplot(data, aes(x=x, y=y)) +
        geom_point() + 
        geom_segment( aes(x=x, xend=x, y=0, yend=y))



p = ggplot(df1,aes(x,y)) + xlab("Sectors by number") + ylab("Job losses") +  geom_segment(data = dfAlElnone,color = "magenta") + geom_point(data = df2,color = "green")  + geom_line(data = dfAlnone,color = "green")  + geom_point(data = dfAlEl.1,color = "red")  + geom_line(data = dfAlEl.1,color = "red") + geom_point(data = dfAlElnone,color = "magenta") + geom_line(data = dfAl.1,color = "blue") +
ggsave("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/jl-full.pdf", height=6, width=8)



Zcomp <- Zbar
pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/heatAcomp.pdf", height=6, width=8)
h <- heatmap.io(Zcomp, RS_label, sectors_x = c(5:34),
                sectors_y = c(7:34), FUN = log, max = 7)
plot(h)
dev.off()


Zcomp <- Zbar
pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/heatfullAcomp.pdf", height=6, width=8)
hfull <- heatmap.io(Zcomp, RS_label, sectors_x = c(5:34,71),
                    sectors_y = c(7:34,71), FUN = log, max = 7)
plot(hfull)
dev.off()


ZexAl <- Zbar
ZexAl[21,] <- rep(0,71)
ZexAl[,21] <- rep(0,71)
pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/heatAexAl.pdf", height=6, width=8)
h <- heatmap.io(ZexAl, RS_label, sectors_x = c(5:34,71),
                sectors_y = c(7:34,71), FUN = log, max = 7)
plot(h)
dev.off()



pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/jl-full.pdf", height=6, width=8)
pAl.1 <- ggplot(data.frame(sojlAl.1), type='h', col='blue', xlab = "Sector, by number", ylab="No. Workers")
pAlEl.1 <- ggplot(data.frame(sojlAlEl.1), type='h', col='red', xlab = "Sector, by number", ylab="No. Workers")
grid.arrange(grobs= list(pAl.1,pAlEl.1),ncol=2)
dev.off()

pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/jl-full.pdf", height=6, width=8)
pAl <- plot(sojlAl.1, col = 'blue',type='h', main='', xlab = "Sector, by number", ylab="No. Workers", height=6, width=8)
pAlEl <- plot(sojlAlEl.1, col = 'red',type='h', main='', xlab = "Sector, by number", ylab="No. Workers", height=6, width=8)
dev.off()
grid.arrange(grobs = c(pAl, pAlEl), ncol = 2)



               dev.off()  
  # re-define data and overwrite top layer inheritance


sum(wcomp)/sum(xcomp)

print(pAlfull <- plot(data.frame(sojlAlfull), type='h', col='blue', main='', xlab = "Sector, by number", ylab="No. Workers"))
print(pAlElfull <- plot(data.frame(sojlAlElfull), type='h', col='red'))
dev.off()
print(p <- plot_grid(pAlfull,pAlElfull,ncol=2))
pdf("~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/jl-full.pdf", height=6, width=8)




plot.histogram

barplot(prop.table(table(sojlAlfull,sojlAlElfull)),beside=T)



library(ggplot2)

sojlAlElnone



.95*Zval[8,21] 
- Eval[21] 
+ .8*Zval[26,21]
(Eval[22]+Eval[23]+Eval[24])
