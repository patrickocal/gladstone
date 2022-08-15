#plot gross value added for each sector
sum(gvacomp- gvaexAlfull)
plot(gvacomp - gvaexAlfull)
plot(gvacomp - gvaexAlnone)

sogvaAlfull <- gvacomp- gvaexAlfull
sogvaAlfull <- sogvaAlfull[-c(73)]
sogvaAlfull <- sogvaAlfull[-c(72)]
sogvaAlfull <- sogvaAlfull[-c(71)]
sogvaAlfull

sogvaAlpar <- gvacomp- gvaexAlpar
sogvaAlpar <- sogvaAlpar[-c(73)]
sogvaAlpar <- sogvaAlpar[-c(72)]
sogvaAlpar <- sogvaAlpar[-c(71)]
sogvaAlpar

sogvaAlnone <- gvacomp- gvaexAlnone
sogvaAlnone <- sogvaAlnone[-c(73)]
sogvaAlnone <- sogvaAlnone[-c(72)]
sogvaAlnone <- sogvaAlnone[-c(71)]
sogvaAlnone

sogvaAlElfull <- gvacomp- gvaexAlElfull
sogvaAlElfull <- sogvaAlElfull[-c(73)]
sogvaAlElfull <- sogvaAlElfull[-c(72)]
sogvaAlElfull <- sogvaAlElfull[-c(71)]
sogvaAlElfull

sogvaAlElpar <- gvacomp- gvaexAlElpar
sogvaAlElpar <- sogvaAlElpar[-c(73)]
sogvaAlElpar <- sogvaAlElpar[-c(72)]
sogvaAlElpar <- sogvaAlElpar[-c(71)]
sogvaAlElpar

sogvaAlElnone <- gvacomp- gvaexAlElnone
sogvaAlElnone <- sogvaAlElnone[-c(73)]
sogvaAlElnone <- sogvaAlElnone[-c(72)]
sogvaAlElnone <- sogvaAlElnone[-c(71)]
sogvaAlElnone

df5 <- data.frame(Scenario=rep(c("BSL alone closes","BSL and GPS close"), each=70),
                  Sectors=rep(1:70,2),
                  Decrease_in_GVA=c(sogvaAlpar,sogvaAlElpar))
df1 <- data.frame(Sectors=rep(1:70,1),
                  Decrease_in_GVA=c(sogvaAlpar))
print(df1)
df4 <- data.frame(Sectors=rep(1:70,1),
                  Decrease_in_GVA=c(sogvaAlElpar))
print(df4)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scale_colour_manual(values=cbbPalette)
dev.off()
p = ggplot(df5,aes(Sectors,Decrease_in_GVA),xlim=.5:70.5) +
  labs(y= "Decrease in Gross Value Added ($M)", x = "Sectors") +
  ggtitle("Decrease in Gross Value Added per Sector for the Central Scenarios") +
  theme(plot.title = element_text(hjust = 0.5,size=22)) +
  geom_segment(aes(x=Sectors,xend=Sectors,y=0,yend=Decrease_in_GVA),color="black",size=1/10,alpha=1/10) +
  geom_point(aes(shape=Scenario),size=1.5)+
  geom_point(aes(color=Scenario),size=.6)+
  geom_line(data = df4, color = "#00C5CD",size=7/10)+
  geom_line(data = df1, color = "#FF8000",size=7/10)+
  scale_colour_hue(l=60)+
  scale_fill_manual(values=cbbPalette)+
  scale_x_continuous(limits = c(.85,70.15), expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 10) +
  guides(color = guide_legend(nrow=1),
         shape = guide_legend(nrow=1))
#  theme_bw()
#  theme(panel.grid = element_blank(),
#        panel.border = element_blank())
plot(p)
ggpar(p,
      #x.text.angle = 60
      font.family="Times",legend="top",font.legend=14,legend.box="vertical")


