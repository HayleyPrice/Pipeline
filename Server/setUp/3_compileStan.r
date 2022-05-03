path = '/mnt/hc-storage/users/hprice/Pipeline/QPROTmodel/QPROTmodel.stan'   ##path to Stan Model


library(rstan)
rstan_options(auto_write = TRUE)
sm <- stan_model(file = path, 
                 save_dso = TRUE)

save('sm', file = 'QPROTmodel.RData')


# /home/hprice/R/x86_64-pc-linux-gnu-library/4.1'



