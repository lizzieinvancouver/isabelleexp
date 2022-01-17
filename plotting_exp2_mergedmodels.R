## Started 16 January 2022 ##
## By Lizzie ##

## Plotting for Stan models for Isabelle's experiments ##
## Combining EACH species level model into one plot ##

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

# Start here and think of sensible way to grab and merge posteriors
# So that code is all ready for PLOTTING ... 
acerpost <- as.matrix(m2.stan.ace.rs)

# I think we need these from each model ... 
acerpost[,"(Intercept)"]
acerpost[,"as.factor(treatment.adj)2_14/26"]
acerpost[,"as.factor(treatment.adj)3_14/26d"]
acerpost[,"as.factor(treatment.adj)4_26/14"]
acerpost[,"as.factor(treatment.adj)5_14/22"]
acerpost[,"as.factor(treatment.adj)6_10/26"]
