
#dataSet <- args[1]
#comp <- args[2]

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTmodel/LocalScripts')

#outFile <- paste(dataSet, "_QLresult_", comp, ".csv", sep = "")
outFile <- "PXD004682_QPROTmodel_testResult.csv"

subsets <- list.files(pattern = "out_*", full.names = FALSE)

results <- NULL

for (subset in subsets) {
  
  data <- tryCatch(read.csv(file = subset), error = function(e) NULL)
  results <- rbind(results, data)
}

# Set the spike in and background proteins
# NB - substring to be searched for in ProteinName
spike <- "HUMAN"
background <- "YEAS"

results$Zstatistic <- abs(results$Zstatistic)
results <- results[order(-results$Zstatistic), ]

results$spike <- mapply(grepl, pattern = spike, x = results$Protein)
results$background <- mapply(grepl, pattern = background, x = results$Protein)
cols <- sapply(results, is.logical)
results[,cols] <- lapply(results[,cols], as.numeric)

results$spike <- mapply(grepl, pattern = spike, x = results$Protein)
results$background <- mapply(grepl, pattern = background, x = results$Protein)
cols <- sapply(results, is.logical)
results[,cols] <- lapply(results[,cols], as.numeric)

# Create running tally
results$spikeTot <- cumsum(results$spike)
results$bgTot <- cumsum(results$background)

# Calculates FDR and sensitivity
results$FDR <- (results$bgTot/(results$spikeTot + results$bgTot))
totSpike1 <- sum(results$spike)
totBG1 <- sum(results$bgTot)
results$Sensitivity <- (results$spikeTot/totSpike1)
results$Sensitivity[is.na(results$Sensitivity)] <- 0
results$Qval <- 0

#results <- results[,2:10]
write.csv(results, file = outFile)


