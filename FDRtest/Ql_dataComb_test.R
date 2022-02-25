library(mclust)
library(factoextra)
library(MASS)
library(ggplot2)
library(FDRestimation)
library(fdrtool)
library(qvalue)
library(locfdr)

#dataSet <- args[1]
#comp <- args[2]

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/FDRtest') 
inFile <- "QPROT_FDRtest"
outFile <- "QPROT_FDRtoolOut.csv"

example_qprot_fdr <- read.table('example_qprot_fdr', header = TRUE, sep = '\t')
example_qprot_fdr <- example_qprot_fdr[order(example_qprot_fdr$Zstatistic), ]
example_qprot_density <- read.table('example_qprot_density', header = TRUE, sep = '\t')

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
#qprot_result <- results

qprot_result <- read.table(inFile, header = TRUE, sep = '\t')
qprot_result <- qprot_result[, c(1:19)]
qprot_result <- qprot_result[order(-abs(qprot_result$Zstatistic)), ]
#qprot_result$fdr <- locfdr(qprot_result$Zstatistic)$fdr
qprot_result$lfdr <- fdrtool(qprot_result$Zstatistic)$lfdr
write.table(qprot_result, 'PXD004682_locFDR.txt', sep = ",")


zeds <- qprot_result
pees <- 2*pnorm(-abs(zeds$Zstatistic))

set.seed(1977)
B = 1000
pStat0 <- matrix(NA, nrow = dim(zeds)[1], ncol = B)





est <- p.fdr(pvalues = zeds$pees, zvalues = zeds$Zstatistic, estim.method = 'storey')
zeds$fdr <- est$fdrs


mc <- densityMclust(qprot_result$Zstatistic, G = 2)
plot(mc, what = 'density', data = qprot_result$Zstatistic)
mclass <- Mclust(qprot_result$Zstatistic, G = 2)
plot(mclass, what = 'classification')
summary(mc, classification = TRUE, parameters = TRUE)
mc$parameters

class <- mc$classification
dens <- mc$density
    
den <- cbind(d1, d2, qprot_result$Zstatistic)



p <- ggplot(den) +
    geom_density(aes(x = fdr), colour = 'red') +
    geom_density(aes(x = fdr_est_2), colour = 'pink') +
    geom_density(aes(x = fdr_est), colour = 'blue') 
    
p

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/PathwayAnalysis')
write.table(qprot_result, outFile, sep = ",")



