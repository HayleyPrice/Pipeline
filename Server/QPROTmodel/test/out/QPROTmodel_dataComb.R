#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

dataSet <- args[1]
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/subs', sep = ""))

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/QPROTmodel/test/out')
norms <- c("None", "Loess-G","RLR-G","VSN-G","TI-G","MedI-G","AI-G","Quantile")

for (norm in norms) {
    
    ##Combine subsets
    subsets <- list.files(pattern = norm, full.names = FALSE)
    
    results <- NULL
    
    for (subset in subsets) {
        
        data <- tryCatch(read.csv(file = subset, header = FALSE), error = function(e) NULL)
        cols <- data[1,]
        data <- data[-1,]
        results <- rbind(results, data)
    }
    
    ## calculate adj p values
    results <- rbind(cols, results)
    #results <- results[order(-abs(results$Zstatistic)), ]
    #results$fdr <- locfdr(results$Zstatistic)$fdr
    #results$pVal <- 2*pnorm(-abs(results$Zstatistic))
    #results$BHpVal <- p.adjust(results$pVal, method = "BH", n = length(results$pVal))
    
    ## print combined reults file for PA analysis
    setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/Results', sep = ""))
    fileOut <- paste(dataset, "_", norm, "_DEresults.csv", sep = '')
    write.csv(results, fileOut, row.names = FALSE)
}



