#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTmodel/ServerScripts')
#args <- c('PXD004682_QPROTin_testShort', 'PXD004682', '5')

#test if there is at least one argument: if not, return an error
if (length(args)==0) {
 stop('At least one argument must be supplied (input file).n', call.=FALSE)
}
dataSet <- args[2]

fileIn <- args[1]
#fileOut <- paste(subset) 

model_input <- read.table(fileIn, sep = '\t', header = TRUE)
model_input<- model_input[, 1:14]

proteins <- nrow(model_input)
split <- as.numeric(args[3])
subs <- ceiling(proteins/split)

#Replaces zero values
model_input[model_input == 0] <- 0.00001
#Logs the intensities
#model_input[,2:7] <- log(model_input[2:7])

setwd('/mnt/hc-storage/users/hprice/Pipeline/QPROTmodel/test')
#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/QPROTmodel/test')
for (i in 1:split) {
 
  firstLine <- (subs * (i-1)) + 1
  #print(firstLine)
  lastLine <- subs * i
  #print(lastLine)
  dataSub <- model_input[firstLine:lastLine,]
  loop <- i
  outFile <- paste(fileIn, '_subset', loop, '.csv', sep = '')
  write.csv(dataSub, outFile)
  
}

setwd('/mnt/hc-storage/users/hprice/Pipeline/QPROTmodel')
