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

df5 <- data.frame(Scenario=rep(c("1a", "2a","1b","2b","1c","2c"), each=70),
                  Sectors=rep(1:70,6),
                  Gross_Value_Added=c(sogvaAlfull,sogvaAlElfull,sogvaAlpar,sogvaAlElpar,sogvaAlnone,sogvaAlElnone))
df1 <- data.frame(Sectors=rep(1:70,1),
                  Gross_Value_Added=c(sogvaAlfull))
print(df1)
df4 <- data.frame(Sectors=rep(1:70,1),
                  Gross_Value_Added=c(sogvaAlElnone))
print(df4)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scale_colour_manual(values=cbbPalette)
dev.off()
p = ggplot(df5,aes(Sectors,Gross_Value_Added),xlim=.5:70.5) +
  geom_segment(aes(x=Sectors,xend=Sectors,y=0,yend=Gross_Value_Added),color="black",size=1/10,alpha=1/10) +
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


