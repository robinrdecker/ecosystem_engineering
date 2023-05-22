
#Support functions for ecosystem engineering

################################################################################
#PLOTTING FUNCTIONS
################################################################################

makePlot <- function(data, byGeneration = 10, title=""){
  #Set color palette
  #colors = matlab.like(max(data$Atoms$Timestep)+1)
  colors = viridis(max(data$Atoms$Timestep)+1)
  #Make initial plot, set up axes
  plotdata <- data$df[ which(data$df$Timestep == 0), ]
  plot(plotdata$Quality, plotdata$Unoccupied, type= "l", xlab = "habitat quality", 
       ylab = "unoccupied area", main = title,
       xlim=c(min(data$df$Quality), max(data$df$Quality)),
       ylim=c(min(min(data$df$Unoccupied), min(data$Atoms$Unoccupied)), max(max(data$df$Unoccupied), max(data$Atoms$Unoccupied))),
       col = colors[1], lwd=2.0)
  atomdata <- data$Atoms[ which(data$Atoms$Timestep == 0), ]
  points(atomdata$Quality, atomdata$Unoccupied, col = colors[1], lwd=2.0)
  
  #Adds all the other distributions, one generation at a time
  for (i in 1:(as.integer(max(data$df$Timestep)/byGeneration))){
    plotdata <- data$df[ which(data$df$Timestep == i*byGeneration), ]
    lines(plotdata$Quality, plotdata$Unoccupied, col = colors[i*byGeneration +1], lwd=2.0)
    atomdata <- data$Atoms[ which(data$Atoms$Timestep == i*byGeneration), ]
    points(atomdata$Quality, atomdata$Unoccupied, col = colors[i*byGeneration +1], lwd=2.0)
  }
}

plotTimestep <- function(data, tmstp = 0, title=""){
  #Make initial plot, set up axes
  plotdata <- data$df[ which(data$df$Timestep == tmstp), ]
  plot(plotdata$Quality, plotdata$Unoccupied, type= "l", xlab = "Quality", 
       ylab = "Area Unoccupied", main = title, 
       xlim=c(min(data$df$Quality), max(data$df$Quality)),
       ylim=c(0, 4))
  atomdata <- data$Atoms[ which(data$Atoms$Timestep == tmstp), ]
  points(atomdata$Quality, atomdata$Unoccupied)
}

plotTimestepAppendix <- function(data, tmstp = 0, title=""){
  #Make initial plot, set up axes
  plotdata <- data$df[ which(data$df$Timestep == tmstp), ]
  plot(plotdata$Quality, plotdata$Unoccupied, type= "l", xlab = "elevation", 
       ylab = "occupied area", main = title, 
       xlim=c(min(data$df$Quality), max(data$df$Quality)),
       ylim=c(0, 4))
  atomdata <- data$Atoms[ which(data$Atoms$Timestep == tmstp), ]
  #points(atomdata$Quality, atomdataOccupied)
  points(-2, 1.5)
  points(1, 2)
}

makePlotO <- function(data, byGeneration = 10, title=""){
  #Set color palette
  #colors = matlab.like(max(data$Atoms$Timestep)+1)
  colors = viridis(max(data$Atoms$Timestep)+1)
  #Make initial plot, set up axes
  plotdata <- data$df[ which(data$df$Timestep == 0), ]
  atomdata <- data$Atoms[ which(data$Atoms$Timestep == 0), ]
  plot(plotdata$Quality, plotdata$Occupied, type= "l", xlab = "habitat quality", 
       ylab = "occupied area", main = title,
       xlim=c(min(data$df$Quality), max(data$df$Quality)),
       ylim=c(min(min(data$Atoms$Occupied),min(data$df$Occupied)), max(max(data$Atoms$Occupied),max(data$df$Occupied))),
       col = colors[1], lwd=2.0)
  points(atomdata$Quality, atomdata$Occupied, col = colors[1], lwd=2.0)
  
  #Adds all the other distributions, one generation at a time
  for (i in 1:(as.integer(max(data$df$Timestep)/byGeneration))){
    plotdata <- data$df[ which(data$df$Timestep == i*byGeneration), ]
    lines(plotdata$Quality, plotdata$Occupied, col = colors[i*byGeneration +1], lwd=2.0)
    atomdata <- data$Atoms[ which(data$Atoms$Timestep == i*byGeneration), ]
    points(atomdata$Quality, atomdata$Occupied, col = colors[i*byGeneration +1], lwd=2.0)
  }
}

plotTimestepO <- function(data, tmstp = 0, title=""){
  #Make initial plot, set up axes
  plotdata <- data$df[ which(data$df$Timestep == tmstp), ]
  plot(plotdata$Quality, plotdata$Occupied, type= "l", xlab = "Quality", 
       ylab = "Area", main = title,
       xlim=c(min(data$df$Quality), max(data$df$Quality)),
       ylim=c(min(data$df$Occupied), max(data$df$Occupied)))
  atomdata <- data$Atoms[ which(data$Atoms$Timestep == tmstp), ]
  points(atomdata$Quality, atomdata$Occupied)
}

################################################################################
#PARALLELIZATION FUNCTIONS
################################################################################

iterate <- function(subintervals, timesteps, d, m, l, mu, sigma,
                    rMax, b, c, initialHabitatQuality,
                    doEngineer, doErosion, doRiseSeaLevel, doGrow){
  
  print(d)
  print(m)
  modelData <- runModelForOccupied(subintervals, 
                                   timesteps, 
                                   d,
                                   m,
                                   l,
                                   mu,
                                   sigma,
                                   rMax,
                                   b,
                                   c,
                                   initialHabitatQuality,
                                   doEngineer,
                                   doErosion,
                                   doRiseSeaLevel,
                                   doGrow) 
  return (modelData)
}

###################################################################
