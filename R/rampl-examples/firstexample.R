firstexample <- function(solver=NULL, modelDirectory=NULL) {
  library(rAMPL)
  # Create an AMPL instance
  ampl <- new(AMPL)

  ## If the AMPL installation directory is not in the system search path:
  # env <- new(Environment, "full path to the AMPL installation directory")
  # ampl <- new(AMPL, env)

  if (!is.null(solver)) {
    ampl$setOption("solver", solver)
  }

  # Read the model and data files.
  if (is.null(modelDirectory)) {
    modelDirectory <- "./models"
  }
  ampl$read(paste(modelDirectory, "/diet/diet.mod", sep=""))
  ampl$readData(paste(modelDirectory, "/diet/diet.dat", sep=""))

  # Solve
  ampl$solve()

  # Get objective entity by AMPL name
  totalcost <- ampl$getObjective("Total_Cost")
  # Print it
  cat(sprintf("objective is: %f\n", totalcost$value()))

  # Reassign data - specific instances
  cost <- ampl$getParameter("cost")
  indices <- c("BEEF", "HAM")
  values <- c(5.01, 4.55)
  cost$setValues(data.frame(indices, values))
  cat(sprintf("Increased costs of beef and ham.\n"))

  # Resolve and display objective
  ampl$solve()
  cat(sprintf("New objective value: %f\n", totalcost$value()))

  # Reassign data - all instances
  elements <- c(3, 5, 5, 6, 1, 2, 5.01, 4.55)
  cost$setValues(elements)
  cat(sprintf("Updated all costs.\n"))

  # Resolve and display objective
  ampl$solve()
  cat(sprintf("New objective value: %f\n", totalcost$value()))

  # Get the values of the variable Buy in a dataframe object
  buy <- ampl$getVariable("Buy")
  df <- buy$getValues()
  # Print them
  print(df)

  # Get the values of an expression into a DataFrame object
  df2 <- ampl$getData("{j in FOOD} 100*Buy[j]/Buy[j].ub")
  # Print them
  print(df2)
}
