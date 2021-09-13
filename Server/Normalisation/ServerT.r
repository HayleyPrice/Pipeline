
# #!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop('At least one argument must be supplied (input file).n', call.=FALSE)
}

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/Normalisation')

fileName <- args[1]
#fileName <- 'subset1'
#fileIn <- paste(fileName, '.csv', sep = '')
outFile <- paste(fileName, "_Tresults.csv", sep = '')

input <- read.table(fileName, sep = '\t', header = TRUE)
input <- na.omit(input)

compN <- input[, c(2:14)]
colnames(compN) <- c("N", "N", "N", "N", "N", "N", "N", "T", "T", "T", "T", "T", "T")

#Variable declaration
ttest<- NULL
pVal <- NULL

#Initialization
count <- (1:nrow(compN))

#Performs t.test
for (k in count){
  index1 <- grep ("N", colnames(compN))
  index2 <- grep ("T", colnames(compN))
  #
  # #check if columns exist
  if(length(index1) == 0) next
  if(length(index2) == 0) next
  
  #Welch Two Sample t-test - EXCEL: two-tailed, unequal variance =T.TEST($C2:$E2,$F2:$H2,2,3)
  ttest <- t.test (as.numeric(compN[k,][index1]), na.rm = TRUE, as.numeric(compN[k,][index2]), na.rm = TRUE, paired = FALSE)
  
  pVal <- rbind (pVal, ttest$p.value)
  
}

# Puts row names on pVals
results <- cbind(input, pVal)

setwd('/mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis/T')
write.csv(results, outFile)



