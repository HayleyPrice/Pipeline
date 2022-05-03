#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

path <- args[1]
stanFile = paste(path, '/QPROTmodel/QPROTmodel.stan', sep = "")   ##path to Stan Model


library(rstan)
rstan_options(auto_write = TRUE)
sm <- stan_model(file = stanFile, 
                 save_dso = TRUE)

save('sm', file = 'QPROTmodel.RData')


# /home/hprice/R/x86_64-pc-linux-gnu-library/4.1'



