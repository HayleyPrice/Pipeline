#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
dataSet <- args[1]
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/QPROTresults/', dataSet, sep = ""))

results <- data.frame(matrix(ncol = 8, nrow = 0))
headers <- c('DE', 'Normalisation',"Threshold", "Go_BP", "Go_MF", "Go_CC",
             "Total", 'DEs')
colnames(results) <- headers

#norms <- c("Loess-G","RLR-G","VSN-G","TI-G","MedI-G","AI-G","Quantile")

#for (norm in norms) {

    #subsets <- list.files(pattern = norm, full.names = FALSE)
    #subsets <- list.files(pattern = '_fdr_PAout', full.names = FALSE)
    subsets <- list.files(pattern = '_QMfdr_PAout', full.names = FALSE)
    
    for (subset in subsets) {
        
        data <- tryCatch(read.csv(file = subset), error = function(e) NULL)
        results <- rbind(results, data)
    }
    
    results <- results[order(-results$Total), ]
    #best <- results[1,]
    # 
    #bestResults <- rbind(bestResults, best)
#}

write.csv(results, 'ResultsSummary_QPROT_test.csv')
