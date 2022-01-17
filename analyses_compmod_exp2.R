## Started 24 October 2019 ##
## By Lizzie ##

## Isabelle's experiments ##

## TO DO NEXT ... ##

############################################
## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# library(lme4)
library(rstan)
library(rstanarm)
# library(plyr)
# library(dplyr)

setwd("~/Documents/git/projects/treegarden/misc/isabelle_expe")

## Get the data ##
## This version has night, DIFF etc. temperatures
d <- read.delim("data_expe2_v2.csv", header=TRUE, sep=";")
d$MaxT <- NA
for(i in c(1:nrow(d))){
    d$MaxT[i] <- max(d$DayT[i], d$NightT[i])
    }

acer <- subset(d, Species=="Acer")

acer.null <- stan_lmer(BBdelay ~ (1|Tree_ID/Cutting_code), data=acer, cores=4)
acer.diff <- stan_lmer(BBdelay ~ DIFF + (1|Tree_ID/Cutting_code), data=acer, cores=4)
acer.max <- stan_lmer(BBdelay ~ MaxT + (1|Tree_ID/Cutting_code), data=acer, cores=4)
acer.diffmax <- stan_lmer(BBdelay ~  MaxT + DIFF + (1|Tree_ID/Cutting_code), data=acer, cores=4)
# summary(acermod)

loo.acer.null <- loo(acer.null, save_psis = TRUE, k_threshold = 0.7)
loo.acer.diff <- loo(acer.diff, save_psis = TRUE, k_threshold = 0.7)
loo.acer.max <- loo(acer.max, save_psis = TRUE, k_threshold = 0.7)
loo.acer.diffmax <- loo(acer.diffmax, save_psis = TRUE, k_threshold = 0.7)

# Least negative is the best ... but need to consider SE when comparing models
compare_models(loo.acer.null, loo.acer.diff, loo.acer.max, loo.acer.diffmax) 
