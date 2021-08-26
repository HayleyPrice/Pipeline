library(rstan)
rstan_options(auto_write = TRUE)
sm <- stan_model(file = '/mnt/hc-storage/users/hprice/Rbatch/QPROTmodel/QPROTmodel.stan', 
                 save_dso = TRUE)

save('sm', file = 'QPROTmodel.RData')


# /home/hprice/R/x86_64-pc-linux-gnu-library/4.1'



