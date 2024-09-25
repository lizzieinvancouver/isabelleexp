## Started 16 January 2022 ##
## By Lizzie ##

## Plotting for Stan models for Isabelle's experiments ##
## Combining EACH species level model into one plot ##

## To do ... ##
# Checking the plot is doing what we want ...
# Did a quick spot check of df versus what was plotted, seems good! #
# Also told Isabelle to check. #

################
## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(rstan) # oddly, this seems to have alpha command!

setwd("~/Documents/git/projects/treegarden/misc/isabelle_expe")

# Get the ful model output 
load("output/m2.stan.ace.rs.Rdata")
load("output/m2.stan.bet.rs.Rdata")
load("output/m2.stan.fag.rs.Rdata")
load("output/m2.stan.que.rs.Rdata")

# Some code while thinking of sensible way to grab and merge posteriors
if(FALSE){
acerpost <- as.matrix(m2.stan.ace.rs)

# I think we need these from each model ... 
acerpost[,"(Intercept)"]
acerpost[,"as.factor(treatment.adj)2_14/26"]
acerpost[,"as.factor(treatment.adj)3_14/26d"]
acerpost[,"as.factor(treatment.adj)4_26/14"]
acerpost[,"as.factor(treatment.adj)5_14/22"]
acerpost[,"as.factor(treatment.adj)6_10/26"]
# 2.5 of acerpost intercept is 33.8
quantile(acerpost[,"(Intercept)"], probs = c(0.025, 0.5)) # yay
}
# And in the end we want ...
# mean, 2.5, 25, 50, 75 and 97.5 for... each species and all species combined

# Okay, we we go!
latbi <- c("acer", "betula", "fagus", "quercus", "overall")

treatcolnames <- c("(Intercept)",
                   "as.factor(treatment.adj)2_14/26",
                   "as.factor(treatment.adj)3_14/26d",
                   "as.factor(treatment.adj)4_26/14",
                   "as.factor(treatment.adj)5_14/22",
                   "as.factor(treatment.adj)6_10/26")

postmean <- (as.matrix(m2.stan.ace.rs)+
                  as.matrix(m2.stan.bet.rs)+
                  as.matrix(m2.stan.fag.rs)+
                  as.matrix(m2.stan.que.rs))/4
                              
posteriorz <- list(acerpost=as.matrix(m2.stan.ace.rs),
                   betpost=as.matrix(m2.stan.bet.rs),
                   fagpost=as.matrix(m2.stan.fag.rs),
                   quepost=as.matrix(m2.stan.que.rs),
                   postmean=postmean)

df <- data.frame(treat=NA,
                 species=NA,
                 mean=NA,
                 twopt5per=NA,
                 twenty5per=NA,
                 fiftyer=NA,
                 seventy5per=NA,
                 ninety7pt5per=NA)


for (postnum in seq_along(posteriorz)){
    posthere <- posteriorz[[postnum]]
      for (whichtreat in treatcolnames){
          dfaddtreat <- data.frame(
                    treat=whichtreat,
                    species=latbi[postnum],
                    mean=mean(posthere[,whichtreat]),
                    twopt5per=quantile(posthere[,whichtreat], probs=c(0.025)),
                    twenty5per=quantile(posthere[,whichtreat], probs=c(0.25)),
                    fiftyer=quantile(posthere[,whichtreat], probs=c(0.5)),
                    seventy5per=quantile(posthere[,whichtreat], probs=c(0.75)),
                    ninety7pt5per=quantile(posthere[,whichtreat], probs=c(0.975)))
          df <- rbind(df, dfaddtreat)
    }
                    
}


## Super, now we need to plot this ...

# First, set up a df with the intercepts set to 0
dfzero <- df
for (spname in seq_along(latbi)){
    for(i in c(3:8)){
        dfzero[which(dfzero$treat=="(Intercept)" & dfzero$species==latbi[spname]),i] <-
            df[which(df$treat=="(Intercept)" & df$species==latbi[spname]),i]-df[which(df$treat=="(Intercept)" & df$species=="overall"),3]
        # Note that you can set df$species=="overall" --> df$species==latbi[spname] and ...
        # then you get each point centered on its own model (but everything piles up then on zero)
}
}

dfzero[which(dfzero$treat=="(Intercept)"),]

# Want muploteachspeciesalone_wavgpost_notzero.pdf? (First version that I made for Isabelle)
# Replace dfzero below with df

# Set up plotting
library(RColorBrewer)
my.pal <- rep(brewer.pal(n = 4, name = "Set1"), 4)
my.pal2 <- c(brewer.pal(n = 9, name = "Set1"), "darkred")

# display.brewer.all()
my.pch <- c(15:18)
cexpch <- 1
alphahere = 0.6 # transparency
alphahere.lighter = 0.2

jitterpt <- -0.2 # push tree ID estimates below main estimate


cil2 <- "twopt5per"
ciu2 <- "ninety7pt5per"
cil <- "twenty5per"
ciu <- "seventy5per"
xlim <- c(-20, 17)
legendxpos <- 9
width <- 7
height <- 6

pdf(file.path("figures/manuscriptfigures/muploteachspeciesalone_wavgpost.pdf"),
    width = width, height = height)
par(xpd=FALSE)
par(mar=c(5,7,3,7))
plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=c(0,6),
     xlab="Model estimate delay in budburst (days)", ylab="",
         main="Experiment 2: Model of each species alone, with average shown")
axis(2, at=1:6, labels=rev(c("20/20", "14/26", "14/26d", "26/14", "14/22", "10/26")), las=1)
abline(v=0, lty=2, col="darkgrey")
for(whichsp in c(1:(length(latbi)-1))){
dfhere <- subset(dfzero, species==latbi[whichsp])
# Okay, we go through the 6 levels of treatment  
for(i in 1:6){ 
  pos.y <- (6:1)[i]
  lines(c(dfhere[i,cil2], 
          dfhere[i,ciu2]), 
          rep(pos.y+jitterpt,2),
          col=alpha(my.pal2[whichsp], alphahere.lighter), type="l", lwd=2)
  lines(c(dfhere[i,cil], 
          dfhere[i,ciu]), 
          rep(pos.y+jitterpt,2),
          col=alpha(my.pal2[whichsp], alphahere), type="l", lwd=2)
  }
pos.y <- (6:1)+jitterpt
pos.x <- dfhere[,"mean"]+jitterpt
points(pos.x, pos.y, cex=1, pch=19, col=alpha(my.pal2[whichsp], alphahere))
}

dfoverall <- subset(dfzero, species=="overall")
 
for(i in 1:6){
   pos.y <- (6:1)
   lines(c(dfoverall[i,cil2], dfoverall[i,ciu2]), rep(pos.y[i], 2),
       col=alpha("black", alphahere.lighter), type="l", lwd=2)
}
    
points(dfoverall[1:6,"mean"], (6:1), cex=cexpch, pch=19, col="black")

for(i in 1:6){
   pos.y <- (6:1)
   lines(c(dfoverall[i,cil], dfoverall[i,ciu]), rep(pos.y[i], 2), col="black",
       type="l", lwd=2)
}

par(xpd=TRUE) # so I can plot legend outside
legend(legendxpos, 1.5, c("Acer", "Betula", "Fagus", "Quercus"),
   pch=19, col=alpha(my.pal2[1:4], alphahere),
   cex=cexpch, bty="n")
dev.off()

write.table(df, file="output/m2.stan_eachspeciesalone_wavgpost.csv", sep=";", row.names=TRUE)

## Update in September 2024
## Isabelle would also like the values for each tree 

# Okay, so here's the parameters 
colMeans(posteriorz[[1]])
# I discussed with Isabelle and she will just take the full set of means (all parameters)
# A bit of a hack, I build the dataframe off the first species and loop through remaining 3
allmodoutputmeanz <- as.data.frame(colMeans(posteriorz[[1]]))
names(allmodoutputmeanz) <- c("mean")
allmodoutputmeanz$species <- latbi[1]
restofthespp <- c()
for(i in c(1:3)){
    nextspp <- as.data.frame(colMeans(posteriorz[[i+1]]))
    names(nextspp) <- c("mean")
    nextspp$species<- latbi[i+1]
    allmodoutputmeanz <- rbind(allmodoutputmeanz, nextspp)
}

write.table(allmodoutputmeanz, file="output/m2.stan_eachspeciesalone_meanoutputallparams.csv", sep=";", 
    row.names=TRUE)

