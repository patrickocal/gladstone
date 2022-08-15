#for full comparison we put all the percentage changes together
#and create tables for latex##
#first output
#gather the percentage changes for each scenario
library(scales)
library(xtable)
p <- c(perchxexAlfull,
  perchxexAlpar,
  perchxexAlnone,
  perchxexAlElfull,
  perchxexAlElpar,
  perchxexAlElnone)
#place into an appropriate matrix
perchxex <- t(matrix(p,3,2))
perchxex
#create a new version of the same matrix, but with % symbols
p <- perchxex/100
p <- percent(p, accuracy = .1, scale = 100, prefix = "",
        suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE)
percchxex <- matrix(p,2,3)
percchxex
#construct table with absolute figures and percentages
#first check the sum of total output is correct
sum(xcomp[1:70])
chxex <- perchxex*sum(xcomp[1:70])/100
chxext <- formatC(perchxex*sum(xcomp[1:70])/100,1,format="f")
comchxex <- matrix(as.character(c(chxext,percchxex)),2,6)
comchxex
#now create labels for the table
clabelcom <- matrix(as.character(c(" (a) Full","(b) Part","(c) None","Full","Part","None")))
rlabel <- matrix(as.character(c("(1) BSL", "(2) BSL & GPS")))
#rename the rows and columns
rownames(comchxex) <- rlabel
colnames(comchxex) <- clabelcom
comchxex
#create a table for latex
res.table <-xtable(comchxex,caption = "Decrease in Total Output (\\$\\textsc{m}) per scenario",
                   align = "lrrr|rrr",digits=c(1,1,1,1,1,1,1))
#view it
print(res.table, scalebox=1, caption.placement = "top",format.args = list(big.mark = ",", decimal.mark = "."))
# output to the home folder for the latex document
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comchxex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))

###now we do the same for wages###

#gather the percentage changes for each scenario
p <- c(perchwexAlfull,
       perchwexAlpar,
       perchwexAlnone,
       perchwexAlElfull,
       perchwexAlElpar,
       perchwexAlElnone)
#place into an appropriate matrix
perchwex <- t(matrix(p,3,2))
perchwex
#create a new version of the same matrix, but with % symbols
p <- perchwex/100
p <- percent(p, accuracy = .1, scale = 100, prefix = "",
             suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE)
percchwex <- matrix(p,2,3)
percchwex
#construct table with absolute figures and percentages
#first check the sum of total employee compensation is correct
sum(wcomp)
chwex <- perchwex*sum(wcomp[1:70])/100
chwext <- formatC(perchwex*sum(wcomp[1:70])/100,1,format="f")
chwex
comchwex <- matrix(as.character(c(chwext,percchwex)),2,6)
comchwex
#now create labels for the table
clabelcom <- matrix(as.character(c(" (a) Full","(b) Part","(c) None","Full","Part","None")))
rlabel <- matrix(as.character(c("(1) BSL", "(2) BSL & GPS")))
#rename the rows and columns
rownames(comchwex) <- rlabel
colnames(comchwex) <- clabelcom
comchwex
#create a table for latex
res.table <-xtable(comchwex,caption = "Decrease in Total Employee compensation (\\$\\textsc{m}) per scenario",
                   align = "lrrr|rrr",digits=1)
#view it
print(res.table, scalebox=1, caption.placement = "top",format.args = list(big.mark = ",", decimal.mark = "."))
# output to the home folder for the latex document
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comchwex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))

#for indirect wages just the matrix
#gather the percentage changes for each scenario
p <- c(indperchwexAlfull,
       indperchwexAlpar,
       indperchwexAlnone,
       indperchwexAlElfull,
       indperchwexAlElpar,
       indperchwexAlElnone)
#place into an appropriate matrix
indperchwex <- t(matrix(p,3,2))
indperchwex

chinwex <- indperchwex*xcomp[71]/100
chinwex

#some calculations
sum(Acomp[,21]) + Mval[21]/xcomp[21]-174.6/xcomp[21]
sum(Acomp[,8]) + Mval[8]/xcomp[8]+733/xcomp[8]

#compute premium Aluminium sector salaries
Alsal <- Zval[71,21]/920
Alsal
xcomp[71]/Alsal
#compute mean salary 
avsal <- Zval[71,21]/920-.01
avsal
lf <- xcomp[71]/avsal
lf
#some checks vs population
sum(gvacomp)/62000
#electricity sector
1/3*Zval[71,26]/Alsal
Zval[71,21]/Alsal

#compute direct job losses from difference of total and indirect change in wages
chdwex <- chwex - chinwex
chdwex

jldex <- chdwex/Alsal
jldex
jldex[2,] <- jldex[2,]-11
jldex
#construct table of percentage changes
perjldex <- jldex/lf*100
perjldex
#add % symbol
p <- perjldex/100
p <- percent(p, accuracy = .1, scale = 100, prefix = "",
             suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE)
percjldex <- matrix(p,2,3)
percjldex
#round up to 0dp and format as text 
jldext <- formatC(jldex,0,format="f")
#construct combined table
comjldex <- matrix(c(jldext,percjldex),2,6)
comjldex
rownames(comjldex) <- rlabel
colnames(comjldex) <- clabelcom
comjldex
#prepare for latex
library(xtable)
res.table <-xtable(comjldex,type="latex",caption = "\\label{comjldex} Direct decrease in employment (\\textsc{fte}) per scenario",
                   align = "lrrr|rrr",digits=1)
res.table
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comjldex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))

#indirect job losses
jlinex <- chinwex/avsal
jlinex
#construct table of percentage changes
perjlinex <- jlinex/lf*100
perjlinex
#add % symbol
p <- perjlinex/100
p <- percent(p, accuracy = .1, scale = 100, prefix = "",
             suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE)
percjlinex <- matrix(p,2,3)
percjlinex
#round up to 1dp
jlinext <- formatC(jlinex,0,format="f")
#construct combined table
comjlinex <- matrix(c(jlinext,percjlinex),2,6)
comjlinex
rownames(comjlinex) <- rlabel
colnames(comjlinex) <- clabelcom
comjlinex
#prepare for latex
library(xtable)
res.table <-xtable(comjlinex,type="latex",caption = "\\label{comjlinex} Indirect decrease in employment (\\textsc{fte})  per scenario",
                   align = "lrrr|rrr",digits=1)
res.table
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comjlinex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))


jltotex <- jldex + jlinex
jltotex
#construct table of percentage changes
perjltotex <- jltotex/lf*100
perjltotex
#add % symbol
p <- perjltotex/100
p <- percent(p, accuracy = .1, scale = 100, prefix = "",
             suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE)
percjltotex <- matrix(p,2,3)
percjltotex
#round up to 1dp
jltotext <- formatC(jltotex,0,format="f")
#construct combined table
comjltotex <- matrix(c(jltotext,percjltotex),2,6)
comjltotex
rownames(comjltotex) <- rlabel
colnames(comjltotex) <- clabelcom
comjltotex
#prepare for latex
library(xtable)
res.table <-xtable(comjltotex,type="latex",caption = "\\label{comjltotex}Total decrease in employment (\\textsc{fte})  per scenario",
                   align = "lrrr|rrr",digits=1)
res.table
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comjltotex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))

jltotex/jldex
#ratio of total to direct job losses
multchwex <- chwex/chdwex
multchwex
#prepare for latex
clabel <- matrix(as.character(c("Full","Part","None")))
rownames(multchwex) <- rlabel
colnames(multchwex) <- clabel
multchwex
library(xtable)
res.table <-xtable(multchwex,type="latex",caption = "\\label{multchwex}Total Employment (\\textsc{fte})  multiplier per scenario",
                   align = "lrrr",digits=1)
res.table
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/multchwex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))

#construct total job losses
#using an average salary of (in millions) 

#construct job losses
#chlex <- round(chwex/avsal,0)
#chlex
chlex <- round(chwex/Alsal,0)
chlex
#check by multiplying by labourforce
#broad
labforbroad <- xcomp[71]/avsal
chlfex <- round(perchwex*labforbroad/100)
chlfex


#construct combined table
comchlex <- matrix(c(chlex,perchwex),2,6)
comchlex
formatC(comchlex, format="d", big.mark=",")
 clabelcom <- matrix(as.character(c(" (a) Full","(b) Part","(c) None","Full","Part","None")))
rlabel <- matrix(as.character(c("(1) BSL", "(2) BSL & GPS")))
rownames(comchlex) <- rlabel
colnames(comchlex) <- clabelcom
comchwex
library(xtable)
res.table <-xtable(comchlex,type="latex",caption = "Total Decrease in Total \\textsc{fte} Employees per scenario",
                   align = "lrrr|rrr",digits=1)
res.table
print(res.table, scalebox=1, caption.placement = "top",format.args = list(big.mark = ",", decimal.mark = "."))
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comchlex.tex", 
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))


###now we do the same for indirect job losses###

#for full comparison we put all the percentage changes together
indperchwex <- t(matrix(c(indperchwexAlfull,
                          indperchwexAlpar,
                          indperchwexAlnone,
                          indperchwexAlElfull,
                          indperchwexAlElpar,
                          indperchwexAlElnone),3,2))
indperchwex
#now create labels for the table
clabel <- matrix(as.character(c("Full","Part","None")))
rlabel <- matrix(as.character(c("(1) BSL", "(2) BSL & GPS")))
#rename the rows and columns
rownames(indperchwex) <- rlabel
colnames(indperchwex) <- clabel
round(indperchwex,1)
#produce a table for latex
library(xtable)
res.table <-xtable(indperchwex,caption = "Percentage change in Total Output for each scenario",
                   align = "rrrr",digits=1)
print(res.table, scalebox=1, caption.placement = "top")

#construct combined table
chinwex <- indperchwex*xcomp[71]/100
chinwext <- formatC(perchwex*sum(wcomp[1:70])/100,1,format="f")
chinwex

comchinwex <- matrix(c(chinwext,indperchwex),2,6)
comchinwex
formatC(comchinwex, format="d", big.mark=",")
 clabelcom <- matrix(as.character(c(" (a) Full","(b) Part","(c) None","Full","Part","None")))
rlabel <- matrix(as.character(c("(1) BSL", "(2) BSL & GPS")))
rownames(comchinwex) <- rlabel
colnames(comchinwex) <- clabelcom
comchinwex
library(xtable)
res.table <-xtable(comchinwex,caption = "Indirect change in Total Employee Compensation per scenario",
                   align = "lrrr|rrr",digits=1)
print(res.table, scalebox=1, caption.placement = "top",format.args = list(big.mark = ",", decimal.mark = "."))
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comchinwex.tex", 
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))






#next gross value added
#gather the percentage changes for each scenario
p <- c(perchgvaexAlfull,
       perchgvaexAlpar,
       perchgvaexAlnone,
       perchgvaexAlElfull,
       perchgvaexAlElpar,
       perchgvaexAlElnone)
#place into an appropriate matrix
perchgvaex <- t(matrix(p,3,2))
perchgvaex
#create a new version of the same matrix, but with % symbols
p <- perchgvaex/100
p <- percent(p, accuracy = .1, scale = 100, prefix = "",
             suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE)
percchgvaex <- matrix(p,2,3)
percchgvaex
#construct table with absolute figures and percentages
#first check the sum of total output is correct
sum(gvacomp)
chgvaex <- round(perchgvaex*sum(gvacomp)/100,1)
chgvaext <- formatC(perchgvaex*sum(gvacomp)/100,1,format="f")
chgvaex
comchgvaex <- matrix(as.character(c(chgvaext,percchgvaex)),2,6)
comchgvaex
#now create labels for the table
clabelcom <- matrix(as.character(c(" (a) Full","(b) Part","(c) None","Full","Part","None")))
rlabel <- matrix(as.character(c("(1) BSL", "(2) BSL & GPS")))
#rename the rows and columns
rownames(comchgvaex) <- rlabel
colnames(comchgvaex) <- clabelcom
comchgvaex
#create a table for latex
res.table <-xtable(comchgvaex,caption = "Total Decrease in Gross Value Added (\\$\\textsc{m}) per scenario",
                   align = "lrrr|rrr",digits=1)
#view it
print(res.table, scalebox=1, caption.placement = "top",format.args = list(big.mark = ",", decimal.mark = "."))
# output to the home folder for the latex document
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comchgvaex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))



#direct effects for gva
chdgvaex <- matrix(c(chdgvaexAlfull,chdgvaexAlElfull,chdgvaexAlpar,chdgvaexAlElpar,chdgvaexAlnone,chdgvaexAlElnone),2,3)
clabel <- matrix(as.character(c("Full","Part","None")))
rownames(chdgvaex) <- rlabel
colnames(chdgvaex) <- clabel
chdgvaex
library(xtable)
res.table <-xtable(chdgvaex,type="latex",caption = "\\label{tab-chdgva} Direct Decrease in Gross Value Added (\\$\\textsc{m}) per scenario",
                   align = "lrrr",digits=1)
res.table
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/chdgvaex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))

#gva multiplier
multgvaex <- matrix(c(multgvaAlfull,multgvaAlElfull,multgvaAlpar,multgvaAlElpar,multgvaAlnone,multgvaAlElnone),2,3)
multgvaex
#prepare for latex
clabel <- matrix(as.character(c("Full","Part","None")))
rownames(multgvaex) <- rlabel
colnames(multgvaex) <- clabel
multgvaex
library(xtable)
res.table <-xtable(multgvaex,type="latex",caption = "\\label{tab-multgva}Gross Value Added multiplier per scenario",
                   align = "lrrr",digits=1)
res.table
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/multgvaex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))




#next exports
#gather the percentage changes for each scenario
p <- c(percheexAlfull,
       percheexAlpar,
       percheexAlnone,
       percheexAlElfull,
       percheexAlElpar,
       percheexAlElnone)
#place into an appropriate matrix
percheex <- t(matrix(p,3,2))
percheex
#create a new version of the same matrix, but with % symbols
p <- percheex/100
p <- percent(p, accuracy = .1, scale = 100, prefix = "",
             suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE)
perccheex <- matrix(p,2,3)
perccheex
#construct table with absolute figures and percentages
#first check the sum of total output is correct
sum(ecomp)
cheex <- round(percheex*sum(ecomp)/100,1)
cheext <- formatC(percheex*sum(ecomp)/100,1,format="f")
cheex
comcheex <- matrix(as.character(c(cheext,perccheex)),2,6)
comcheex
#now create labels for the table
clabelcom <- matrix(as.character(c(" (a) Full","(b) Part","(c) None","Full","Part","None")))
rlabel <- matrix(as.character(c("(1) BSL", "(2) BSL & GPS")))
#rename the rows and columns
rownames(comcheex) <- rlabel
colnames(comcheex) <- clabelcom
comcheex
#create a table for latex
res.table <-xtable(comcheex,caption = "Decrease in Total Exports (\\$\\textsc{m}) per scenario",
                   align = "lrrr|rrr",digits=1)
#view it
print(res.table, scalebox=1, caption.placement = "top",format.args = list(big.mark = ",", decimal.mark = "."))
# output to the home folder for the latex document
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comcheex.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))

lf
jltotex*26000/lf

chgvaex*26000/lf
chgvaex
gladfactor <- 27000/lf
sum(xcomp[1:70])*gladfactor
sum(gvacomp)*gladfactor


chelgva <- gvacomp[26]-matrix (c(gvaexAlfull[26],gvaexAlElfull[26],gvaexAlpar[26],gvaexAlElpar[26],gvaexAlnone[26],gvaexAlElnone[26]),2,3)
perchelgva <- chelgva/gvacomp[26]
p <-perchelgva
p <- percent(p, accuracy = .1, scale = 100, prefix = "",
             suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE)
p <- matrix(p,2,3)
p
chelgvat <- formatC(chelgva,1,format="f")
chelgvat
comchelgvat <- matrix(as.character(c(chelgvat,p)),2,6)
comchelgvat
#now create labels for the table
clabelcom <- matrix(as.character(c(" (a) Full","(b) Part","(c) None","Full","Part","None")))
rlabel <- matrix(as.character(c("(1) BSL", "(2) BSL & GPS")))
#rename the rows and columns
rownames(comchelgvat) <- rlabel
colnames(comchelgvat) <- clabelcom
comchelgvat
#create a table for latex
res.table <-xtable(comchelgvat,caption = "Total Decrease in Electricity Sector Value Added (\\$\\textsc{m}) per scenario",
                   align = "lrrr|rrr",digits=1)
#view it
print(res.table, scalebox=1, caption.placement = "top",format.args = list(big.mark = ",", decimal.mark = "."))
# output to the home folder for the latex document
print.xtable(res.table,file = "~/OneDrive - The University of Queensland/_QTC/_Gladstone/_tex/comchelgvat.tex",
             caption.placement = "top",
             table.placement = "H",
             format.args = list(big.mark = ",", decimal.mark = "."))


