#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop('At least one argument must be supplied (input file).n', call.=FALSE)
}

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/Normalisation')
norm <- args[1]
dataSet <- args[2]

fileName <- paste(dataSet, "_" , norm, "-normalized_Log2", sep = "")
#fileName <- 'subset1'
#fileIn <- paste(fileName, '.csv', sep = '')
outFile <- paste(fileName, "_Tresults.csv", sep = '')

input <- as.matrix(read.table(fileName, sep = '\t', header = FALSE))
#input <- na.omit(input)

cols <- input[1,]
input <- input[-1,]
colNum <- length(cols)

compN <- input[, c(2:colNum)]
colnms <- cols[-1]
#colnms <- as.vector(as.character(colnames))
colnames(compN) <- colnms

cols <- c(cols, 'pVal', 'BHpVal')

#Variable declaration
ttest<- NULL
pVal <- NULL

#Initialization
count <- (1:nrow(compN))

#Performs t.test
for (k in count){
  index1 <- grep ("1", colnames(compN))
  index2 <- grep ("2", colnames(compN))
  #
  # #check if columns exist
  if(length(index1) == 0) next
  if(length(index2) == 0) next
  
  #Welch Two Sample t-test - EXCEL: two-tailed, unequal variance =T.TEST($C2:$E2,$F2:$H2,2,3)
  ttest <- t.test (as.numeric(compN[k,][index1]), na.rm = TRUE, as.numeric(compN[k,][index2]), na.rm = TRUE, paired = FALSE)
  
  pVal <- rbind (pVal, ttest$p.value)
  
}

#pVal <- rbind('pVal', pVal)
# Puts row names on pVals
results <- cbind(input, pVal)
#results <- results[order(pVal), ]
BHpVal <- p.adjust(as.numeric(pVal), method = "BH", n = length(pVal))
results <- cbind(input, pVal, BHpVal)
results <- results[order(pVal), ]
colnames(results) <- cols
setwd('/mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis/T')

write.csv(results, outFile)



