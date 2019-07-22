## Started 15 July 2019 ##
## By Lizzie ##

## Plotting for Stan models for Isabelle's experiments ##

################
## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(rstan) # oddly, this seems to have alpha command!

setwd("~/Documents/git/projects/treegarden/isabelle_expe")

# Get the models (probably need to change to tab delim!)
# read.delim("data_expe1.csv", header=TRUE, sep=";")
m2.stan.ace.rs.df <- read.table("output/m2.stan.ace.rs.df.csv", sep=";")
m2.stan.bet.rs.df <- read.delim("output/m2.stan.bet.rs.df.csv", sep=";")
m2.stan.fag.rs.df <- read.delim("output/m2.stan.fag.rs.df.csv", sep=";")
m2.stan.que.rs.df <- read.delim("output/m2.stan.que.rs.df.csv", sep=";")
nameshere <- c("mean", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "mcse")
names(m2.stan.ace.rs.df) <- nameshere
names(m2.stan.bet.rs.df) <- nameshere
names(m2.stan.fag.rs.df) <- nameshere
names(m2.stan.que.rs.df) <- nameshere

## Look at one model
m2.stan.ace.rs.df[1:10,]

### Set up plotting
library(RColorBrewer)
my.pal <- rep(brewer.pal(n = 4, name = "Set1"), 4)
my.pal2 <- c(brewer.pal(n = 9, name = "Set1"), "darkred")

# display.brewer.all()
my.pch <- c(15:18)
alphahere = 0.6 # transparency
alphahere.lighter = 0.2

jitterpt <- -0.2 # push tree ID estimates below main estimate


#################################
## f(x)s to make plotting easier 

# function to plot for each model with tree ID including actual intercept value
plot.randslopes <- function(model, spname, nameforfig, height, width, cil, ciu, xlim, legendxpos){
    modelhere <- model
pdf(file.path(paste("figures/muplot", spname, nameforfig, ".pdf", sep="")),
    width = width, height = height)
par(xpd=FALSE)
par(mar=c(5,7,3,7))
plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=c(0,6),
     xlab=paste("Model estimate delay in budburst (", spname, ")", sep=""), ylab="",
         main=paste("Experiment 2: Model for", spname))
axis(2, at=1:6, labels=rev(c("20/20 (intcpt)", "14/26", "14/26d", "26/14", "14/22", "10/26")), las=1)
abline(v=0, lty=2, col="darkgrey")
treelist <- c("Tree:1", "Tree:2", "Tree:3", "Tree:4","Tree:5",
    "Tree:6", "Tree:7", "Tree:8", "Tree:9", "Tree:10")
for(whichtree in c(1:length(treelist))){
dfhere <- modelhere[grep(whichtree[1], rownames(modelhere)),][1:6,]
for(i in 1:6){
  pos.y <- (6:1)[i]
  lines(c(dfhere[i,cil]+modelhere[i,cil], dfhere[i,ciu]+modelhere[i,ciu]), rep(pos.y+jitterpt,2), col=alpha(my.pal2[whichtree], alphahere), type="l", lwd=2)
 # lines(c(dfhere[i,"2.5%"], dfhere[i,"97.5%"]), rep(pos.y,2), col=alpha(my.pal2[whichtree], alphahere.lighter), type="l", lwd=2)
  }
pos.y <- (6:1)+jitterpt
pos.x <- dfhere[,"mean"]+modelhere[1:6,"mean"]+jitterpt
points(pos.x, pos.y, cex=1, pch=19, col=alpha(my.pal2[whichtree],alphahere))
}

points(modelhere[1:6,"mean"], (6:1), cex=1, pch=19, col="black")
for(i in 1:6){
   pos.y <- (6:1)
   lines(c(modelhere[i,cil], modelhere[i,ciu]), rep(pos.y[i], 2), col="black",
       type="l", lwd=2)
}

par(xpd=TRUE) # so I can plot legend outside
legend(legendxpos, 6, c("Tree:1", "Tree:2", "Tree:3", "Tree:4","Tree:5",
    "Tree:6", "Tree:7", "Tree:8", "Tree:9", "Tree:10"),
   pch=19, col=alpha(my.pal2[1:10], alphahere),
   cex=0.75, bty="n")
dev.off()
}


# function to plot for each model with tree ID with intercept set at zero
if(FALSE){ # Below I used for debugging code, but also explains what the function needs to run!
model <- m2.stan.ace.rs.df # the dataframe of the model 
spname <- "Acer" # what species? Shows up in legends
nameforfig <- "whatever" # part of your file name
height <- 7 # height of the figure
width <- 6 # width of the figue
cil <- "2.5%" # credible interval lower
ciu <- "97.5%" # credible interval upper
xlim <- c(-20, 20) # the x axis limits
legendxpos <- 23 # where to put the legend (must be al little bigger than max of xlim) 
} 
plot.randslopes.zeroint <- function(model, spname, nameforfig, height, width, cil, ciu, xlim, legendxpos){
    modelhere <- model
    modelhere[1,1:6] <- modelhere[1,1:6]-modelhere[1,1]
pdf(file.path(paste("figures/muplot", spname, nameforfig, ".pdf", sep="")),
    width = width, height = height)
par(xpd=FALSE)
par(mar=c(5,7,3,7))
plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=c(0,6),
     xlab=paste("Model estimate delay in budburst (", spname, ")", sep=""), ylab="",
         main=paste("Experiment 2: Model for", spname))
axis(2, at=1:6, labels=rev(c("20/20 (intcpt)", "14/26", "14/26d", "26/14", "14/22", "10/26")), las=1)
abline(v=0, lty=2, col="darkgrey")
treelist <- c("Tree:1", "Tree:2", "Tree:3", "Tree:4","Tree:5",
    "Tree:6", "Tree:7", "Tree:8", "Tree:9", "Tree:10")
for(whichtree in c(1:length(treelist))){
dfhere <- modelhere[grep(whichtree[1], rownames(modelhere)),][1:6,]
for(i in 1:6){
  pos.y <- (6:1)[i]
  lines(c(dfhere[i,cil]+modelhere[i,cil], dfhere[i,ciu]+modelhere[i,ciu]), rep(pos.y+jitterpt,2), col=alpha(my.pal2[whichtree], alphahere), type="l", lwd=2)
 # lines(c(dfhere[i,"2.5%"], dfhere[i,"97.5%"]), rep(pos.y,2), col=alpha(my.pal2[whichtree], alphahere.lighter), type="l", lwd=2)
  }
pos.y <- (6:1)+jitterpt
pos.x <- dfhere[,"mean"]+modelhere[1:6,"mean"]+jitterpt
points(pos.x, pos.y, cex=1, pch=19, col=alpha(my.pal2[whichtree],alphahere))
}
points(modelhere[1:6,"mean"], (6:1), cex=1, pch=19, col="black")
for(i in 1:6){
   pos.y <- (6:1)
   lines(c(modelhere[i,cil], modelhere[i,ciu]), rep(pos.y[i], 2), col="black",
       type="l", lwd=2)
}

par(xpd=TRUE) # so I can plot legend outside
legend(legendxpos, 6, c("Tree:1", "Tree:2", "Tree:3", "Tree:4","Tree:5",
    "Tree:6", "Tree:7", "Tree:8", "Tree:9", "Tree:10"),
   pch=19, col=alpha(my.pal2[1:10], alphahere),
   cex=0.75, bty="n")
dev.off()
}


#################################
## Make the plots!

# old versions
plot.randslopes(m2.stan.ace.rs.df, "Acer", "exp2randslopes.50per", 7, 6, "25%", "75%", c(-10, 46), 48)
plot.randslopes(m2.stan.bet.rs.df, "Betula", "exp2randslopes.50per", 7, 6, "25%", "75%", c(-10, 46), 48)
plot.randslopes(m2.stan.fag.rs.df, "Fagus", "exp2randslopes.50per", 7, 6, "25%", "75%", c(-20, 56), 60)
plot.randslopes(m2.stan.que.rs.df, "Quercus", "exp2randslopes.50per", 7, 6, "25%", "75%", c(-10, 46), 48)

# intercept centered at zero
plot.randslopes.zeroint(m2.stan.ace.rs.df, "Acer", "exp2randslopes.noint.50per", 7, 6, 
    "25%", "75%", c(-20, 20), 23)
plot.randslopes.zeroint(m2.stan.ace.rs.df, "Acer", "exp2randslopes.noint.95per", 7, 6,
    "2.5%", "97.5%", c(-20, 20), 23)
plot.randslopes.zeroint(m2.stan.bet.rs.df, "Betula", "exp2randslopes.noint.50per", 7, 6,
    "25%", "75%", c(-20, 20), 23)
plot.randslopes.zeroint(m2.stan.bet.rs.df, "Betula", "exp2randslopes.noint.95per", 7, 6,
    "2.5%", "97.5%", c(-20, 20), 23)
plot.randslopes.zeroint(m2.stan.fag.rs.df, "Fagus", "exp2randslopes.noint.50per", 7, 6,
    "25%", "75%", c(-30, 20), 23)
plot.randslopes.zeroint(m2.stan.fag.rs.df, "Fagus", "exp2randslopes.noint.95per", 7, 6,
    "2.5%", "97.5%", c(-30, 20), 23)
plot.randslopes.zeroint(m2.stan.que.rs.df, "Quercus", "exp2randslopes.noint.50per", 7, 6,
    "25%", "75%", c(-20, 20), 23)
plot.randslopes.zeroint(m2.stan.que.rs.df, "Quercus", "exp2randslopes.noint.95per", 7, 6,
     "2.5%", "97.5%", c(-20, 20), 23)

