#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTmodel/ServerScripts')
#args <- c('PXD004682_QPROTin_testShort', 'PXD004682', '5')

#test if there is at least one argument: if not, return an error
if (length(args)==0) {
 stop('At least one argument must be supplied (input file).n', call.=FALSE)
}
norm <- args[1]
dataSet <- args[2]
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, "/NormalisedData", sep = ""))

fileIn <- paste(dataSet, "_" , norm, "-normalized_Log2", sep = "")
#fileOut <- paste(subset) 
#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/Normalisation')
#fileIn <- 'PXD004682_AI-G-normalized_Log2'
model_input <- read.table(fileIn, sep = '\t', header = FALSE)
cols <- model_input[1,]
colNum <- length(cols)
model_input<- model_input[-1, 1:colNum]

proteins <- nrow(model_input)
split <- as.numeric(args[3])
subs <- ceiling(proteins/split)

#Replaces zero values
model_input[model_input == 0] <- 0.00001
#Logs the intensities
#model_input[,2:7] <- log(model_input[2:7])

setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, "/Qmodel/subs", sep = ""))

for (i in 1:split) {
 
  firstLine <- (subs * (i-1)) + 1
  #print(firstLine)
  lastLine <- subs * i
  #print(lastLine)
  dataSub <- model_input[firstLine:lastLine,]
  loop <- i
  outFile <- paste(fileIn, '_subset', loop, '.csv', sep = '')
  dataSub <- rbind(cols,dataSub)
  write.csv(dataSub, outFile)
  
}

