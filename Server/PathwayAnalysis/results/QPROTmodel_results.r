
#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/PathwayAnalysis/results')

results <- data.frame(matrix(ncol = 9, nrow = 0))
headers <- c('DE', 'Normalisation',"Threshold", "Go_BP", "Go_MF", "Go_CC",
             "Total", 'DEs')
colnames(results) <- headers

norms <- c("Loess-G","RLR-G","VSN-G","TI-G","MedI-G","AI-G","Quantile")

for (norm in norms) {
    
    #subsets <- list.files(pattern = norm, full.names = FALSE)
    #subsets <- list.files(pattern = 'PAout', full.names = FALSE)
    
    #for (subset in subsets) {
        
        data <- tryCatch(read.csv(file = paste("PAout_", norm, sep = "")), error = function(e) NULL)
        data <- data[, c(2:9)]
        results <- rbind(results, data)
    #}
    

    #best <- results[1,]
    # 
    #bestResults <- rbind(bestResults, best)
}
results <- results[order(-results$Total), ]
write.csv(results, 'ResultsSummary.csv')
