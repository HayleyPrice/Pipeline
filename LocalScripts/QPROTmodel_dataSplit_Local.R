#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
dataSet <- args[2]

model_input <- read.table(args[1], sep = "\t", header = TRUE)
model_input<- model_input[, 1:14]

proteins <- nrow(model_input)
split <- as.numeric(args[3])
subs <- ceiling(proteins/split)

#Replaces zero values
model_input[model_input == 0] <- 0.00001
#Logs the intensities
#model_input[,2:14] <- log(model_input[2:14])

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTmodel/LocalScripts/test')

for (i in 1:split) {
   
  firstLine <- (subs * (i-1)) + 1
  #print(firstLine)
  lastLine <- subs * i
  #print(lastLine)
  dataSub <- model_input[firstLine:lastLine,]
  loop <- i
  outFile <- paste("subset", loop, ".csv", sep = "")
  write.csv(dataSub, outFile)
  
}

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTmodel/LocalScripts')
