library(rstan)
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
#library(Rcpp)

# #!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

start <- Sys.time()
setwd('/mnt/hc-storage/users/hprice/Rbatch/QPROTmodel/test')

ch = 4
c = 1
loc_input = 650
qprot_iterations = loc_input*c
wu = loc_input/2
AD = 0.9

fileName <- args[1]
#fileName <- "subset1"
fileIn <- paste(fileName, ".csv", sep = "")
ouputFile <- paste("out_", fileName, ".csv", sep = "")

sink(ouputFile, append=TRUE)
#sink(ouputFile, append=TRUE, type="message")

QPROTmodel_input <- read.table(fileName, sep = ",", header = TRUE)
#head(QPROTmodel_input)
QPROTmodel_input <- QPROTmodel_input[, 2:15]

#Replaces zero values
QPROTmodel_input[QPROTmodel_input == 0] <- 0.00001
QPROTmodel_input <- na.omit(QPROTmodel_input)


indicator <- c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)             #Describes experimental set up
numProteins <- nrow(QPROTmodel_input)               #Number of proteins in dat-set
numRuns <- length(indicator)                   #Total number of runs

qlone_output <- QPROTmodel_input
qlone_output$LogFC <- 1
qlone_output$Zstatistic <- 1

# # #Performs stan model on each protein in data-set
for (protein in 1:nrow(QPROTmodel_input)) {
  
  intensities <- NULL
  d <- NULL
  
  p <- QPROTmodel_input[protein, 1]
  
  intensities <- as.numeric(QPROTmodel_input[protein, c(2:14)])
  
  #Input data for Stan model
  protData <- list(
    T = indicator,                             #Identifies condition intesity belongs to (0 or 1)
    N = numRuns,                               #Number of intensities per protein
    y = intensities                                 #Protein intensity values, logged
  )
  
  #Runs stan model from file
  fit <- stan(
    #Path to stan model script
    file = '/mnt/hc-storage/users/hprice/Rbatch/QPROTmodel/QPROTmodel.stan',
    data = protData,
    chains = ch,                                #Number of MCMC chains to perform
    warmup = wu,                     #Number of iterations for burnin
    iter = loc_input,                          #Number of iterations on each chain
    cores = c,
    control = list(adapt_delta = AD,
                   max_treedepth = 15)
  )
  
  #Extracts summary data for DE
  d <- summary(fit, pars = c("d"))
  #print(d)
  
  #Extracts the mean and standard DE
  mean <- d$summary[1]
  sd <- d$summary[3]
  
  #Calcuates the log fold change and the Zstatistic for DE
  qlone_output$LogFC[protein] <- d$summary[1]
  qlone_output$Zstatistic[protein] <- d$summary[1]/d$summary[3]
}

write.csv(qlone_output, outfile)

timeTaken <- Sys.time() - start
print(timeTaken)

# Restore output to console
sink() 
sink(type="message")


