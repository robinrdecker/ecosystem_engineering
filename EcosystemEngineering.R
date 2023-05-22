#Ecosystem Engineering

#This code generates the figures  from the manuscript

rm(list = ls())

library(Rcpp)
library(colorRamps) #For colorful plots
library(parallel)
library(lattice)
library(viridis)
library(gridGraphics)
grab_grob <- function(){
  grid.echo()
  grid.grab()
}


setwd() 

sourceCpp("spartina.cpp")
#The C++ script depends on several proprietary header files created by Numerical
#Recipes 3, which I have not included. If you have access to these, you can 
#include them to compile the program. Otherwise you can use your own favorite
#Romberg integration routine, or any other numerical quadrature you prefer.

source("EcosystemEngineeringFunctions.R")

################################################################################
# EXAMPLES
################################################################################

DataA <- runModel(subintervals = 3000, #Need 2000 for accuracy
                  timesteps = 50, 
                  d = 0.1, #deterioration, erosion
                  m = 0.1, #habitat modification, engineering
                  l = 2.1, #effect of nearby individuals on engineering ability
                  mu = 0.0, #value of optimal habitat quality
                  sigma = 0.5, #std dev about that optimal habitat quality
                  rMax = 2.0, #maximum growth rate (at optimal quality)
                  b = 0.4, #dispersal ability of propagules 0.4
                  c = 0, #rate of sea level rise
                  initialHabitatQuality = 2,
                  doEngineer = 0,
                  doErosion = 1,
                  doRiseSeaLevel =0,
                  doGrow = 0)

DataB <- runModel(subintervals = 3000, #Need 2000 for accuracy
                  timesteps = 50, 
                  d = 0.1, #deterioration, erosion
                  m = 0.1, #habitat modification, engineering
                  l = 2.1, #effect of nearby individuals on engineering ability
                  mu = 0.0, #value of optimal habitat quality
                  sigma = 0.5, #std dev about that optimal habitat quality
                  rMax = 2.0, #maximum growth rate (at optimal quality)
                  b = 0.4, #dispersal ability of propagules 0.4
                  c = 0, #rate of sea level rise
                  initialHabitatQuality = 2,
                  doEngineer = 1,
                  doErosion = 0,
                  doRiseSeaLevel =0,
                  doGrow = 0)

DataC <- runModel(subintervals = 3000, #Need 2000 for accuracy
                  timesteps = 30, 
                  d = 0.0, #deterioration, erosion
                  m = 0.0, #habitat modification, engineering
                  l = 2.1, #effect of nearby individuals on engineering ability
                  mu = 0.0, #value of optimal habitat quality
                  sigma = 0.5, #std dev about that optimal habitat quality
                  rMax = 2.0, #maximum growth rate (at optimal quality)
                  b = 0.4, #dispersal ability of propagules 0.4
                  c = 0.1, #rate of sea level rise
                  initialHabitatQuality = 2,
                  doEngineer = 1,
                  doErosion = 1,
                  doRiseSeaLevel =1,
                  doGrow = 0)

DataD <- runModel(subintervals = 3000, #Need 2000 for accuracy
                  timesteps = 30, 
                  d = 0.1, #deterioration, erosion
                  m = 0.1, #habitat modification, engineering
                  l = 2.1, #effect of nearby individuals on engineering ability
                  mu = 0.0, #value of optimal habitat quality
                  sigma = 0.5, #std dev about that optimal habitat quality
                  rMax = 2.0, #maximum growth rate (at optimal quality)
                  b = 0.4, #dispersal ability of propagules 0.4
                  c = 0, #rate of sea level rise
                  initialHabitatQuality = 2,
                  doEngineer = 0,
                  doErosion = 0,
                  doRiseSeaLevel =0,
                  doGrow = 1)

par(cex.axis=1.5)
par(cex.lab=2.5)
par(oma=c(0,2,0,0)) #outer margins (bottom, left, top, right)
par(mar=c(5.1, 5.0, 4.1, 2.1)) #margins (5.1, 4.1, 4.1, 2.1)

makePlot(DataA, byGeneration = 3, "")
dev.copy(pdf, 'Fig 1A.pdf', width=9, height=6);dev.off()

makePlot(DataB, byGeneration = 3, "")
dev.copy(pdf, 'Fig 1B.pdf', width=9, height=6);dev.off()

makePlot(DataC, byGeneration = 3, "")
dev.copy(pdf, 'Fig 1C.pdf', width=9, height=6);dev.off()

makePlotO(DataD, byGeneration = 3, "")
dev.copy(pdf, 'Fig 1D.pdf', width=9, height=6);dev.off()


plotTimestepO(DataA, tmstp = 0, "Growth")
plotTimestep(DataA, tmstp = 7, "Erosion")

#Appendix
plotTimestepAppendix(DataA, tmstp = 7, "")

################################################################################
#FIGURE 2
################################################################################

#FIGURE 2A

engineerRates <- c(3.0, 6.0, 0.0)


fig1ListLow <- mclapply(engineerRates, iterate, 
                        subintervals = 3000,
                        timesteps = 100,
                        d = 0.2, #deterioration, erosion
                        l = 0.1, #effect of nearby individuals on engineering ability
                        mu = 0.0, #value of optimal habitat quality
                        sigma = 0.3, #std dev about that optimal habitat quality
                        rMax = 3.0, #maximum growth rate (at optimal quality)
                        b = 0.4, #dispersal ability of propagules 0.4
                        c = 0, #rate of sea level rise
                        initialHabitatQuality = 0, #Initially suboptimal
                        doEngineer = 1,
                        doErosion = 1,
                        doRiseSeaLevel =0,
                        doGrow = 1)

fig1MatLow <- matrix(unlist(fig1ListLow), ncol=length(engineerRates), nrow=101)


dev.off()

x11(width=4, height=5)
times <- 0:100
par(bg=NA, las=1) 
plot(times, fig1MatLow[,1], type='l', xlab="", ylab="area occupied",
     ylim=c(0.0,1.0), lwd=1.5, cex.axis=0.8) 
title(xlab = "year", mgp = c(2, 1, 0))
lines(times, fig1MatLow[,2], lty=3, lwd=1.5) 
lines(times, fig1MatLow[,3], lty=5, lwd=1.5) 
legend(60, 0.85,c("m = 0.0", "m = 3.0", "m = 6.0"),
       lty=c(5, 1, 3), cex=0.8, y.intersp=1.5)

plot1 <- grab_grob()

#FIGURE 2C

dvalues <- seq(0, 0.5, 0.025)
mvalues <- seq(0, 10.0, 0.5)

fig2List <- mcmapply(iterate, rep(dvalues,each=length(mvalues)), mvalues,
                     subintervals = 3000,
                     timesteps = 100,
                     l = 0.1, #effect of nearby individuals on engineering ability
                     mu = 0.0, #value of optimal habitat quality
                     sigma = 0.3, #std dev about that optimal habitat quality
                     rMax = 3.0, #maximum growth rate (at optimal quality)
                     b = 0.4, #dispersal ability of propagules 0.4
                     c = 0, #rate of sea level rise
                     initialHabitatQuality = 0, #Initially suboptimal
                     doEngineer = 1,
                     doErosion = 1,
                     doRiseSeaLevel =0,
                     doGrow = 1)

fig2ListLow <- readRDS("fig2List_lowL.rds")
fig2MatLow <- readRDS("fig2Mat_lowL.rds")

xvalues <- mvalues
yvalues <- dvalues

grid <- expand.grid(x=xvalues, y=yvalues)

grid$z <- fig2MatLow[40,] 

plot3 <- levelplot(z~x*y, grid, xlab="engineering (m)", 
                   ylab="erosion (d)", col.regions  = viridis(100))


#FIGURE 2B

engineerRates <- c(0.22, 0.4, 0.0)

fig1ListHigh <- mclapply(engineerRates, iterate,
                         subintervals = 3000,
                         timesteps = 100,
                         d = 0.2, #deterioration, erosion
                         l = 2.1, #effect of nearby individuals on engineering ability
                         mu = 0.0, #value of optimal habitat quality
                         sigma = 0.3, #std dev about that optimal habitat quality
                         rMax = 3.0, #maximum growth rate (at optimal quality)
                         b = 0.4, #dispersal ability of propagules 0.4
                         c = 0, #rate of sea level rise
                         initialHabitatQuality = 0, #Initially suboptimal
                         doEngineer = 1,
                         doErosion = 1,
                         doRiseSeaLevel =0,
                         doGrow = 1)

fig1MatHigh <- matrix(unlist(fig1ListHigh), ncol=length(engineerRates), nrow=101)

dev.off()

x11(width=4, height=5)
times <- 0:100
par(bg=NA, las=1) 
plot(times, fig1MatHigh[,1], type='l', xlab="", ylab="area occupied", 
     lwd=1.5, ylim=c(0.0, 1.0), cex.axis=0.8) 
title(xlab = "year", mgp = c(2, 1, 0))
lines(times, fig1MatHigh[,2], lty=3, lwd=1.5) 
lines(times, fig1MatHigh[,3], lty=5, lwd=1.5) 
legend(60, 0.85,c("m = 0.00", "m = 0.22", "m = 0.40"),
       lty=c(5, 1, 3), cex=0.8, y.intersp=1.5)

plot2 <- grab_grob()

#FIGURE 2D

dvalues <- seq(0, 0.4, 0.025)
mvalues <- seq(0, 1.0, 0.05)

fig2ListHigh <- mcmapply(iterate, rep(dvalues,each=length(mvalues)), mvalues,
                         subintervals = 3000,
                         timesteps = 100,
                         l = 2.1, #effect of nearby individuals on engineering ability
                         mu = 0.0, #value of optimal habitat quality
                         sigma = 0.3, #std dev about that optimal habitat quality
                         rMax = 3.0, #maximum growth rate (at optimal quality)
                         b = 0.4, #dispersal ability of propagules 0.4
                         c = 0, #rate of sea level rise
                         initialHabitatQuality = 0, #Initially suboptimal
                         doEngineer = 1,
                         doErosion = 1,
                         doRiseSeaLevel =0,
                         doGrow = 1)

fig2ListHigh <- readRDS("fig2List_highL.rds")
fig2MatHigh <- readRDS("fig2Mat_highL.rds")

xvalues <- mvalues
yvalues <- dvalues

grid <- expand.grid(x=xvalues, y=yvalues)

grid$z <- fig2MatHigh[40,] 

plot4 <- levelplot(z~x*y, grid, xlab="engineering (m)", 
                   ylab="erosion (d)", col.regions  = viridis(100))

library(gridExtra)
grid.arrange(plot1,plot2,plot3,plot4, ncol=2)
grid.text(c("(a)","(b)", "(c)", "(d)"), x=c(0.05, 0.55, 0.05, 0.55), y=c(0.95, 0.95, 0.47, 0.47), vjust=1, hjust=0, gp=gpar(fontface=1))


################################################################################
#FIGURE 3
################################################################################

# FIGURE 3A

engineerRates <- c(0.22, 0.4, 0.0)

fig3aList <- mclapply(engineerRates, iterate,
                      subintervals = 3000,
                      timesteps = 100,
                      d = 0.2, #deterioration, erosion
                      l = 2.1, #effect of nearby individuals on engineering ability
                      mu = 0.0, #value of optimal habitat quality
                      sigma = 0.3, #std dev about that optimal habitat quality
                      rMax = 3.0, #maximum growth rate (at optimal quality)
                      b = 0.4, #dispersal ability of propagules 0.4
                      c = 0.1, #rate of sea level rise
                      initialHabitatQuality = 0, #Initially suboptimal
                      doEngineer = 1,
                      doErosion = 1,
                      doRiseSeaLevel =1,
                      doGrow = 1)

fig3aMat <- matrix(unlist(fig3aList), ncol=length(engineerRates), nrow=101)

dev.off()

x11(width=4, height=5)
times <- 0:100
par(bg=NA, las=1) 
plot(times, fig3aMat[,1], type='l', xlab="", ylab="area occupied", 
     lwd=1.5, ylim=c(0.0, 1.0), cex.axis=0.8) 
title(xlab = "year", mgp = c(2, 1, 0))
lines(times, fig3aMat[,2], lty=3, lwd=1.5) 
lines(times, fig3aMat[,3], lty=5, lwd=1.5) 
legend(60, 0.85,c("m = 0.00", "m = 0.22", "m = 0.40"),
       lty=c(5, 1, 3), cex=0.8, y.intersp=1.5)

plot3a <- grab_grob()

#FIGURE 3B

cvalues <- seq(0, 0.4, 0.025)
mvalues <- seq(0, 1.0, 0.05)

fig3List <- mcmapply(iterate, rep(mvalues,each=length(cvalues)), cvalues,
                     subintervals = 3000,
                     timesteps = 50,
                     d = 0.2,
                     l = 2.1, #effect of nearby individuals on engineering ability
                     mu = 0.0, #value of optimal habitat quality
                     sigma = 0.3, #std dev about that optimal habitat quality
                     rMax = 3.0, #maximum growth rate (at optimal quality)
                     b = 0.4, #dispersal ability of propagules 0.4
                     initialHabitatQuality = 0, #Initially suboptimal
                     doEngineer = 1,
                     doErosion = 1,
                     doRiseSeaLevel =1,
                     doGrow = 1)

fig3Mat <- matrix(unlist(fig3List), ncol=length(mvalues)*length(cvalues), byrow=FALSE)

xvalues <- cvalues
yvalues <- mvalues

grid <- expand.grid(x=xvalues, y=yvalues)

grid$z <- fig3Mat[40,] 

plot3b <- levelplot(z~x*y, grid, xlab="sea-level rise (c)", 
                    ylab="engineering (m)", col.regions  = viridis(100))

library(gridExtra)
grid.arrange(plot3a,plot3b, ncol=1)
grid.text(c("(a)","(b)"), x=c(0.05, 0.05), y=c(0.95, 0.47), vjust=1, hjust=0, gp=gpar(fontface=1))

################################################################################
# FIGURE 4
################################################################################

#Figure 4A

initialhq <- c(0,1,2)

fig4aList <- mclapply(initialhq, iterate,
                     subintervals = 3000,
                     timesteps = 50, 
                     d = 0.1, #deterioration, erosion
                     m = 0.2,
                     l = 2.1, #effect of nearby individuals on engineering ability
                     mu = 0.0, #value of optimal habitat quality
                     sigma = 0.3, #std dev about that optimal habitat quality
                     rMax = 2.0, #maximum growth rate (at optimal quality)
                     b = 0.4, #dispersal ability of propagules 0.4
                     c = 0, #rate of sea level rise
                     doEngineer = 1,
                     doErosion = 1,
                     doRiseSeaLevel =0,
                     doGrow = 1)

fig4aMat <- matrix(unlist(fig4aList), ncol=3, nrow=51)

dev.off()

x11(width=4, height=5)
times <- 0:50
par(bg=NA, las=1) 
plot(times, fig4aMat[,1], type='l', xlab="", ylab="area occupied", 
     lwd=1.5, ylim=c(0.0, 1.0), cex.axis=0.8) 
title(xlab = "year", mgp = c(2, 1, 0))
lines(times, fig4aMat[,2], lty=3, lwd=1.5) 
lines(times, fig4aMat[,3], lty=5, lwd=1.5) 
legend(20, 0.95,c("suboptimal", "neutral", "optimal"),
       lty=c(1, 3, 5), cex=0.8, y.intersp=1.5)

plot4a <- grab_grob()

#Figure 4B

initialhq <- c(0,1,2)

fig4bList <- mclapply(initialhq, iterate,
                     subintervals = 3000,
                     timesteps = 50, 
                     d = 0.1, #deterioration, erosion
                     m = 0.2,
                     l = 2.1, #effect of nearby individuals on engineering ability
                     mu = 0.0, #value of optimal habitat quality
                     sigma = 0.3, #std dev about that optimal habitat quality
                     rMax = 2.0, #maximum growth rate (at optimal quality)
                     b = 0.4, #dispersal ability of propagules 0.4
                     c = 0.1, #rate of sea level rise
                     doEngineer = 1,
                     doErosion = 1,
                     doRiseSeaLevel =1,
                     doGrow = 1)

fig4bMat <- matrix(unlist(fig4bList), ncol=3, nrow=51)

dev.off()

x11(width=4, height=5)
times <- 0:50
par(bg=NA, las=1) 
plot(times, fig4bMat[,1], type='l', xlab="", ylab="area occupied", 
     lwd=1.5, ylim=c(0.0, 1.0), cex.axis=0.8) 
title(xlab = "year", mgp = c(2, 1, 0))
lines(times, fig4bMat[,2], lty=3, lwd=1.5) 
lines(times, fig4bMat[,3], lty=5, lwd=1.5) 
legend(20, 0.5,c("suboptimal", "neutral", "optimal"),
       lty=c(1, 3, 5), cex=0.8, y.intersp=1.5)

plot4b <- grab_grob()


library(gridExtra)
grid.arrange(plot4a,plot4b, ncol=2)
grid.text(c("(a)","(b)"), x=c(0.05, 0.55), y=c(0.95, 0.95), vjust=1, hjust=0, gp=gpar(fontface=1))


