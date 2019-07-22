## Started 10 July 2018 ##
## By Lizzie ##

## Isabelle's experiments ##

## TO DO NEXT ... ##

############################################
## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## Questions (for both experiments):
# Are there really enough to estimate effects?
# Is flask coded for time in treatment? It should be ... it cannot change for each bud.
# Should treatment be a factor? Yes, I think so (i.e., not divided into day and night temperatures)

library(lme4)
library(rstan)
library(rstanarm)
library(plyr)
library(dplyr)

# Change below if you want to run the Stan models (need multicore)
runstanmodels.e1 <- FALSE
runstanmodels.e2 <- FALSE
runstanmodels.e2.bysp <- TRUE
plotstanmodels.e2.bysp <- TRUE # need to have run the stan models and have the output 

setwd("~/Documents/git/projects/treegarden/isabelle_expe")

## Exp 1 ##
d <- read.delim("data_expe1.csv", header=TRUE, sep=";")

## summarizing data ... the data structure is critical to how lme4 (or related commands)
# will build the model, so it's important to check it seems the data as we want
table(d$species, d$cutting.unq)
table(d$species, d$treatment, d$batch)
table(d$species, d$treatment, d$batch, d$cutting)
table(d$species, d$treatment, d$cutting, d$bud)

# Note however, that flushing was low
dfl <- subset(d, flushed==1)
nrow(dfl)/nrow(d) # 58% flushed
table(dfl$species, dfl$cutting.unq)
# Probably too low to fit cutting

length(unique(d$flask)) # 18 flasks in exp 1?
length(unique(d$cutting))
length(unique(d$bud))

hist(d$BBdelay)

############
## Expe 1 ##
############

# Some models ... 
aov(BBdelay~treatment*species*batch+flask+Error(species/cutting/bud), data=d)
# This is very bad, I would not use AOV for such a complicated model
# "In general, aov() is simpler and faster, but only works for a simpler subset of models." https://stats.stackexchange.com/questions/154700/are-aov-with-error-same-as-lmer-of-lme4-package-in-r

m1 <- lmer(BBdelay~as.factor(treatment)*species*as.factor(batch) +
    (1|species/cutting), data=subset(d, species!="Fagus"))
    # Fagus only in one batch
    # Check: https://stats.stackexchange.com/questions/70556/have-i-correctly-specified-my-model-in-lmer

# m1.old <- lmer(BBdelay~as.factor(treatment)*species*as.factor(batch) +
   # (1|flask)+(1|species/cutting), data=d)
    # Model is still overspecified: "Model is nearly unidentifiable"
    # This is poorly fit, should try MCMC approach

if(runstanmodels.e1){ # run multicore!
m1.stan <- stan_lmer(BBdelay ~ (as.factor(treatment)*species*as.factor(batch))
    + (1|flask) + (1|species/cutting), data = subset(d, species!="Fagus"))
summary(m1.stan)
save(m1.stan, file="m1.exp1.Rdata")
}


# Skip cutting and run non-ME (mixed-effects) model
m1 <- lm(BBdelay~as.factor(treatment)*species*as.factor(batch), data=d)
m1.nfag <- lm(BBdelay~as.factor(treatment)*species*as.factor(batch), data=subset(d, species!="Fagus"))

# summary(m1)

############
## Expe 2 ##
############
d2 <- read.delim("data_expe2.csv", header=TRUE, sep=";")
# d2 <- read.delim("data_expe2_v2.csv", header=TRUE, sep=";") # not using yet ... would need to adjust code a bit, so prefer not to (but may want to merge in for temperature info)

# now get a unique cutting for each tree
d2$cutting.unq <- paste(d2$Flask_code, d2$cutting) 
# check
d2summ <-
      ddply(d2, c("Species", "Tree"), summarise,
      length(unique(cutting.unq)))


## re-sort the treatments
# here's how the model will resort the treatments
sort(unique(d2$Treatment))
d2$treatment.adj <- d2$Treatment
d2$treatment.adj[d2$Treatment=="20/20"] <- "1_20/20"
d2$treatment.adj[d2$Treatment=="14/26"] <- "2_14/26"
d2$treatment.adj[d2$Treatment=="14/26d"] <- "3_14/26d"
d2$treatment.adj[d2$Treatment=="26/14"] <- "4_26/14"
d2$treatment.adj[d2$Treatment=="14/22"] <- "5_14/22"
d2$treatment.adj[d2$Treatment=="oct-26"] <- "6_10/26" # clean-up
sort(unique(d2$treatment.adj))

## 
hist(d2$BBdelay)
hist(log(d2$BBdelay))

## Now run the models ... all species in one model
m2.allspp <- lmer(BBdelay~as.factor(treatment.adj)*Species + (1|Species/Tree), data=d2)
# m2.allspp.n <- lme(BBdelay~as.factor(treatment.adj), random=~1|Species/Tree, data=d2) # same as above buy code for librray(nlme), which gives p-values if you want them

# run with multicore!
if(runstanmodels.e2){
m2.stan <- stan_lmer(BBdelay~as.factor(treatment.adj)*Species +
    (1|Species/Tree), data=d2, cores=4)
m2sum <- summary(m2.stan)
# 647 divergent transitions, maybe because of dup. treatment of species
# so I don't report this model

m2.stan.alt <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (1|Species/Tree), data=d2, cores=4)
m2sum.alt <- summary(m2.stan.alt)
# 4 divergent transitions
m2.rs.stan <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (as.factor(treatment.adj)|Species/Tree), data=d2, cores=4)
# 16 divergent transitions
summary(m2.rs.stan)
}


### Set up plotting
library(RColorBrewer)
my.pal <- rep(brewer.pal(n = 4, name = "Set1"), 4)
my.pal2 <- c(brewer.pal(n = 9, name = "Set1"), "darkred")

# display.brewer.all()
my.pch <- c(15:18)
alphahere = 0.6 # transparency
alphahere.lighter = 0.2

####
if(runstanmodels.e2){
# load("output/m2.stan.alt.Rdata")
# m2sum.alt <- summary(m2.stan.alt)

# ugly, but functional
rownames(m2sum.alt)[1:6] <- paste(rownames(m2sum.alt)[1:6], "Acer")
# loop over species 
par(xpd=FALSE)
par(mar=c(5,7,3,10))
plot(x=NULL,y=NULL, xlim=c(-10, 40), yaxt='n', ylim=c(0,6),
     xlab="Model estimate delay in BB", ylab="", main="Experiment 2")
axis(2, at=1:6, labels=rev(c("20/20", "10/26", "14/22", "14/26", "14/26d", "26/14")), las=1)
abline(v=0, lty=2, col="darkgrey")

modelsplist <- c("Acer", "Betula", "Fagus", "Quercus")
for(whichsp in c(1:4)){
dfhere <- m2sum.alt[grep(modelsplist[whichsp], rownames(m2sum.alt)),][1:6,]
for(i in 1:6){
  pos.y <- (6:1)[i]
  lines(c(dfhere[i,"25%"], dfhere[i,"75%"]), rep(pos.y,2), col=alpha(my.pal[whichsp], alphahere), type="l", lwd=2)
  lines(c(dfhere[i,"2.5%"], dfhere[i,"97.5%"]), rep(pos.y,2), col=alpha(my.pal[whichsp], alphahere.lighter),
     type="l", lwd=2)
  }
pos.y <- (6:1)
pos.x <- dfhere[,"mean"]
points(pos.x, pos.y, cex=1, pch=19, col=alpha(my.pal[whichsp],alphahere))
}

par(xpd=TRUE) # so I can plot legend outside
legend(21, 6, c("Acer", "Betula", "Fagus", "Quercus"), pch=19, col=alpha(my.pal[1:4], alphahere),
   cex=0.75, bty="n")
}

####
## Now run the models ... one species at a time
####
ace <- subset(d2, Species=="Acer")
bet <- subset(d2, Species=="Betula")
fag <- subset(d2, Species=="Fagus")
que <- subset(d2, Species=="Quercus")

# check
table(ace$Treatment, ace$Tree)
table(bet$Treatment, bet$Tree)
table(fag$Treatment, fag$Tree)
table(que$Treatment, que$Tree)

m2.ace <- lmer(BBdelay~as.factor(treatment.adj) + (1|Tree), data=ace)
m2.bet <- lmer(BBdelay~as.factor(treatment.adj) + (1|Tree), data=bet)
m2.fag <- lmer(BBdelay~as.factor(treatment.adj) + (1|Tree), data=fag)
m2.que <- lmer(BBdelay~as.factor(treatment.adj) + (1|Tree), data=que)

# plot just the means of lmer models.... (could add SE)
par(xpd=FALSE)
par(mar=c(5,7,3,7))
plot(x=NULL,y=NULL, xlim=c(-15, 52), yaxt='n', ylim=c(0,6),
     xlab="Model estimate delay in BB", ylab="", main="Experiment 2, separate models")
axis(2, at=1:6, labels=rev(c("20/20", "10/26", "14/22", "14/26", "14/26d", "26/14")), las=1)
abline(v=0, lty=2, col="darkgrey")
pos.y <- (6:1)
pos.xa <- fixef(m2.ace)
pos.xb <- fixef(m2.bet)
pos.xf <- fixef(m2.fag)
pos.xq <- fixef(m2.que)
points(pos.xa, pos.y, cex=1, pch=19, col=alpha(my.pal[1],alphahere))
points(pos.xb, pos.y, cex=1, pch=19, col=alpha(my.pal[2],alphahere))
points(pos.xf, pos.y, cex=1, pch=19, col=alpha(my.pal[3],alphahere))
points(pos.xq, pos.y, cex=1, pch=19, col=alpha(my.pal[4],alphahere))
par(xpd=TRUE) # so I can plot legend outside
legend(21, 6, c("Acer", "Betula", "Fagus", "Quercus"), pch=19, col=alpha(my.pal[1:4], alphahere),
   cex=0.75, bty="n")


if(runstanmodels.e2.bysp){
columnstoextract <- c("mean", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "mcse")
# random intercepts (just for Acer for now)
m2.stan.ace.ri <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (1|Tree), data=ace)
# ... add cutting on random intercepts (all species)
m2.stan.ace.ri.wcut <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (1|Tree/cutting), data=ace)
m2.stan.bet.ri.wcut <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (1|Tree/cutting), data=bet)
m2.stan.fag.ri.wcut <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (1|Tree/cutting), data=fag)
m2.stan.que.ri.wcut <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (1|Tree/cutting), data=que)
# random slopes and intercepts:
# This assumes that you have multiple observations PER tree PER treatment
# 1-2 divergent transitions per model
m2.stan.ace.rs <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (as.factor(treatment.adj)|Tree), data=ace, cores=4)
save(m2.stan.ace.rs, file="output/m2.stan.ace.rs.Rdata")
m2.stan.ace.rs.df <- as.data.frame(summary(m2.stan.ace.rs)[,columnstoextract])
write.table(m2.stan.ace.rs.df, file="output/m2.stan.ace.rs.df.csv", sep=";", row.names=TRUE)

m2.stan.bet.rs <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (as.factor(treatment.adj)|Tree), data=bet, cores=4)
save(m2.stan.bet.rs, file="output/m2.stan.bet.rs.Rdata")
m2.stan.bet.rs.df <- as.data.frame(summary(m2.stan.bet.rs)[,columnstoextract])
write.table(m2.stan.bet.rs.df, file="output/m2.stan.bet.rs.df.csv", sep=";", row.names=TRUE)

m2.stan.fag.rs <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (as.factor(treatment.adj)|Tree), data=fag, cores=4)
save(m2.stan.fag.rs, file="output/m2.stan.fag.rs.Rdata")
m2.stan.fag.rs.df <- as.data.frame(summary(m2.stan.fag.rs)[,columnstoextract])
write.table(m2.stan.fag.rs.df, file="output/m2.stan.fag.rs.df.csv", sep=";", row.names=TRUE)

m2.stan.que.rs <- stan_lmer(BBdelay~as.factor(treatment.adj) +
    (as.factor(treatment.adj)|Tree), data=que, cores=4)
save(m2.stan.que.rs, file="output/m2.stan.que.rs.Rdata")
m2.stan.que.rs.df <- as.data.frame(summary(m2.stan.que.rs)[,columnstoextract])
write.table(m2.stan.que.rs.df, file="output/m2.stan.que.rs.df.csv", sep=";", row.names=TRUE)
    
}


if(plotstanmodels.e2.bysp){
load("output/m2.stan.ace.rs.Rdata")
load("output/m2.stan.bet.rs.Rdata")
load("output/m2.stan.fag.rs.Rdata")
load("output/m2.stan.que.rs.Rdata")

modelhere.full <- m2.stan.ace.rs # cheap! But replace this for the other three species to change plot (watch out for title)
modelhere <- summary(modelhere.full)

# To do ...
# Make intercept 0 
# 20/20
# 14/26
# 14/26d
# 26/14
# 14/22
# 10/26

jitterpt <- -0.2 # push tree ID estimates below main estimate
# loop over trees 
par(xpd=FALSE)
par(mar=c(5,7,3,7))
plot(x=NULL,y=NULL, xlim=c(-10, 46), yaxt='n', ylim=c(0,6),
     xlab="Model estimate delay in BB", ylab="", main="Experiment 2: Model for Acer")
axis(2, at=1:6, labels=rev(c("20/20 (intercept)", "10/26", "14/22", "14/26", "14/26d", "26/14")), las=1)
abline(v=0, lty=2, col="darkgrey")

modelsplist <- c("Tree:1", "Tree:2", "Tree:3", "Tree:4","Tree:5",
    "Tree:6", "Tree:7", "Tree:8", "Tree:9", "Tree:10")
for(whichtree in c(1:length(modelsplist))){
dfhere <- modelhere[grep(modelsplist[whichtree], rownames(modelhere)),][1:6,]
for(i in 1:6){
  pos.y <- (6:1)[i]
  lines(c(dfhere[i,"25%"]+modelhere[i,"25%"], dfhere[i,"75%"]+modelhere[i,"75%"]), rep(pos.y+jitterpt,2), col=alpha(my.pal2[whichtree], alphahere), type="l", lwd=2)
 # lines(c(dfhere[i,"2.5%"], dfhere[i,"97.5%"]), rep(pos.y,2), col=alpha(my.pal2[whichtree], alphahere.lighter), type="l", lwd=2)
  }
pos.y <- (6:1)+jitterpt
pos.x <- dfhere[,"mean"]+modelhere[1:6,"mean"]+jitterpt
points(pos.x, pos.y, cex=1, pch=19, col=alpha(my.pal2[whichtree],alphahere))
}

points(modelhere[1:6,"mean"], (6:1), cex=1, pch=19, col="black")
for(i in 1:6){
   pos.y <- (6:1)
   lines(c(modelhere[i,"25%"], modelhere[i,"75%"]), rep(pos.y[i], 2), col="black",
       type="l", lwd=2)
}

par(xpd=TRUE) # so I can plot legend outside
legend(48, 6, c("Tree:1", "Tree:2", "Tree:3", "Tree:4","Tree:5",
    "Tree:6", "Tree:7", "Tree:8", "Tree:9", "Tree:10"),
   pch=19, col=alpha(my.pal2[1:10], alphahere),
   cex=0.75, bty="n")


}


####
## Old code for recoding by cutting
####
if(FALSE){
# recode cutting, even flasks have cutting 1, odd flasks have cutting 2
d2$flaskID <- as.numeric(unlist(lapply( strsplit(d2$Flask_code, split="-"), "[", 2)))
fleven <- seq(2, 20, by=2)
flodd <- seq(1, 19, by=2)
d2$cutting <- NA
d2$cutting[which(d2$flaskID %in% fleven)] <- 2
d2$cutting[which(d2$flaskID %in% flodd)] <- 1
}
