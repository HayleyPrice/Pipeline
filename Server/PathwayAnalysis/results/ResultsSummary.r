#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

dataset <- args[1]
#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/PathwayAnalysis/results')

results <- data.frame(matrix(ncol = 9, nrow = 0))
headers <- c('Normalisation',"Threshold", "Go_BP", "Go_MF", "Go_CC",
             "Total", 'DEs', 'DE')
colnames(results) <- headers

norms <- c("Loess-G","RLR-G","VSN-G","TI-G","MedI-G","AI-G","Quantile")

for (norm in norms) {
        setwd('/mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis/results')
        data <- tryCatch(read.csv(file = paste(dataset, "_PAout_", norm, ".csv", sep = "")), error = function(e) NULL)
        data <- data[, c(2:9)]
        data$DE <- 'QM'
        results <- rbind(results, data)
}
for (norm in norms) {
    setwd('/mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis/T/results')
    data <- tryCatch(read.csv(file = paste(dataset, "_T_PAout_", norm, ".csv", sep = "")), error = function(e) NULL)
    data <- data[, c(2:9)]
    data$DE <- 'T'
    results <- rbind(results, data)
}

results <- results[order(-results$Total), ]
setwd('/mnt/hc-storage/users/hprice/Pipeline')
write.csv(results, paste(dataset, '_ResultsSummary.csv', sep = ''))
