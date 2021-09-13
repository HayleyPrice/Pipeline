library(ashr)
library(ggplot2)

#dataSet <- args[1]
#comp <- args[2]

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTmodel/LocalScripts/test') 

outFile <- "PXD004682_QPROTmodel_testResult_ENNIX.csv"

# subsets <- list.files(pattern = "out_*", full.names = FALSE)
# 
# results <- NULL
# 
# for (subset in subsets) {
#   
#   data <- tryCatch(read.csv(file = subset), error = function(e) NULL)
#   results <- rbind(results, data)
# }
# 
# zeds <- results
# zeds$pees <- 2*pnorm(-abs(zeds$Zstatistic))
# 
# est <- p.fdr(pvalues = zeds$pees, zvalues = zeds$Zstatistic, estim.method = 'storey')
# zeds$fdr <- est$fdrs

qprot_result <- read.table('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTout_test', header = TRUE, sep = '\t')






p <- ggplot(qprot_result) +
    geom_density(aes(x = fdr), colour = 'red') +
    geom_density(aes(x = fdr_est_2), colour = 'pink') +
    geom_density(aes(x = fdr_est), colour = 'blue') 
    
p

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/PathwayAnalysis')
write.table(qprot_result, "QPROT_FDRtest", sep = "\t")



