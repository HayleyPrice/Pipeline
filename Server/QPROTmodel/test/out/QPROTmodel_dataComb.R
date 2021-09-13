#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/QPROTmodel/test/out')
norms <- c("None", "Loess-G","RLR-G","VSN-G","TI-G","MedI-G","AI-G","Quantile")

for (norm in norms) {
    
    ##Combine subsets
    setwd('/mnt/hc-storage/users/hprice/Pipeline/QPROTmodel/test/out')
    subsets <- list.files(pattern = norm, full.names = FALSE)
    
    results <- NULL
    
    for (subset in subsets) {
        
        data <- tryCatch(read.csv(file = subset), error = function(e) NULL)
        results <- rbind(results, data)
    }
    
    ## calculate adj p values
    results <- results[order(-abs(results$Zstatistic)), ]
    results$pVal <- 2*pnorm(-abs(results$Zstatistic))
    results$BHpVal <- p.adjust(results$pVal, method = "BH", n = length(results$pVal))
    
    ## print combined reults file for PA analysis
    setwd('/mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis')
    fileOut <- paste(norm, "_DEresults.csv", sep = '')
    write.csv(results, fileOut)
}

