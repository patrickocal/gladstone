rm(list = ls())
library(rAMPL)
# Create an AMPL instance
ampl <- new(AMPL)

# path to model file.
modeldir <- "../ampl/code/"
# Read the model file.
ampl$read(paste(modeldir, "maiwar.mod", sep=""))
# path to data
setparval <- function(param_name, datadir="../julia/data/"){
  amplparentity <- ampl$getParameter(param_name)
  parvaldf <- amplparentity$getValues()
  numrow <- nrow(parvaldf)
  newdatadf <- read.csv(file = paste(param_name, ".csv", sep=""))
  newdatadf[newdatadf<0] <- 0
  newdatacol <- matrix(t(as.matrix(newdatadf[,-1])),
                       nrow=numrow, ncol=1,
                       byrow=TRUE)
  parvaldf[, param_name] <- newdatacol
  amplparentity$setValues(parvaldf)
}

setparval("RAW_INV_FLW")
stop
ampl$solve()
# Print out the result
cat(sprintf("Objective: %f\n", ampl$getObjective("pres_disc_val")$value()))
# Set solver, otherwise will default to knitro
#ampl$setOption("solver", solver)


#phi <- ampl$getSet("Sectors")

#sectorNames <- phi$getValues()

#sectorNames$Sectors[9] <- "dI"

#phi$setValues(sectorNames)

#phi$getValues()

#raw_inv_flw <- ampl$getParameter("RAW_INV_FLW")

#raw_inv_flw_val <- raw_inv_flw$getValues()

#raw_con_flw <- ampl$getParameter("RAW_CON_FLW")

#raw_con_flw_val <- raw_con_flw$getValues()

#raw_inv_flw_val[raw_inv_flw_val=="dI"] <- "I"

#data_frame <- raw_inv_flw_val

#df <-data_frame[order(data_frame$index0, data_frame$index1, data_frame$index2),]

#flowsData <- read.csv(file = "capitalFlowsRAS.csv")

#flowsData[flowsData<0] <- 0

#mdat <- matrix(t(as.matrix(flowsData[,2:21])), nrow=400, ncol=1, byrow=TRUE)

#raw_inv_flw_val$RAW_INV_FLW <- mdat

#raw_inv_flw$setValues(raw_inv_flw_val)

#raw_inv_flw$getValues()

#cost$setValues(data.frame(indices, values))

# Just to test that a basic solve as is works:

# Read in a .csv file of data. note that ideally in making this data we want 
# each parameter set to be in one column, of length of the set it runs over, for 
# instance PHI_ADJ would be one column in the csv, of length numSec
#sampleData <- read.csv(file = "sampleData.csv", header=FALSE)

# Pick out single parameter set
#paramData1 <- sampleData$V1

# Set sectors (below is just an example, obviously does not have all the 
# sectors we will be using)
#sectorNames <- c("A", "B", "C", "D", "E", "F", "G", "H")

# Set the values for the set Sectors, and for the parameter PHI_ADJ
# For each parameter you set, you want to include the set that it is indexed 
# over, as below
#ampl$setData(data.frame(Sectors=sectorNames, PHI_ADJ=paramData1), 1, "Sectors")

# I'm leaving this in from dietmodel.R in the examples but commented out, 
# as its a good example of how to import a set indexed over multiple sets

# We would want to pull this in from a .csv and convert to a matrix rather than 
# making it manually as below
# amounts = rbind(
#   c( 60,    8,   8,  40,   15,  70,   25,   60),
#   c( 20,    0,  10,  40,   35,  30,   50,   20),
#   c( 10,   20,  15,  35,   15,  15,   25,   15),
#   c( 15,   20,  10,  10,   15,  15,   15,   10),
#   c(928, 2180, 945, 278, 1182, 896, 1329, 1397),
#   c(295,  770, 440, 430,  315, 400,  379,  450)
# )

# As above we need to include the set that the data is indexed accros
#dimnames(amounts) <- list(nutrients, foods)

# Once again, the matrix has to be converted into a data.frame before we can 
# pass it like this
# df <- data.frame(as.table(amounts))
# colnames(df) <- c("NUTR", "FOOD", "amt")
# # Set the values for the parameter "amt"
# ampl$setData(df, 2, "")

# Solve the model
#ampl$solve()



# Example of how to retreive variables, retained from dietmodel.R 
# in the examples
# df <- ampl$getVariable("Buy")$getValues()
#print(df)

