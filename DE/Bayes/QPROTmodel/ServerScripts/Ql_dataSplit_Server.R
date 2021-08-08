#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#test if there is at least one argument: if not, return an error
if (length(args)==0) {
 stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
dataSet <- args[2]
comp <- args[3]

wd <- paste(dataSet, "/", comp, sep = "")
setwd(wd)

qlone_input <- read.table(args[1], sep = "\t", header = TRUE)
qlone_input <- qlone_input[, 1:7]

proteins <- nrow(qlone_input)
split <- as.numeric(args[4])
subs <- ceiling(proteins/split)

#Replaces zero values
qlone_input[qlone_input == 0] <- 0.01
#Logs the intensities
qlone_input[,2:7] <- log(qlone_input[2:7])

for (i in 1:split) {
   
  firstLine <- (subs * (i-1)) + 1
  #print(firstLine)
  lastLine <- subs * i
  #print(lastLine)
  dataSub <- qlone_input[firstLine:lastLine,]
  loop <- i
  outFile <- paste("subsets/subset", loop, ".csv", sep = "")
  write.csv(dataSub, outFile)
  
}

