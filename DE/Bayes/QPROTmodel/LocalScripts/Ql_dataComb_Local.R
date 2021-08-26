library(mclust)
library(factoextra)
library(MASS)
library(ggplot2)
#dataSet <- args[1]
#comp <- args[2]

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTmodel/LocalScripts/test') 

outFile <- "PXD004682_QPROTmodel_testResult.csv"

subsets <- list.files(pattern = "out_*", full.names = FALSE)

results <- NULL

for (subset in subsets) {
  
  data <- tryCatch(read.csv(file = subset), error = function(e) NULL)
  results <- rbind(results, data)
}

zeds <- results[,c(2,17)]

dens <- densityMclust(zeds$Zstatistic)
plot(dens, what = 'density' )

mc <- Mclust(zeds$Zstatistic, G = 2)
plot(mc, what = 'classification')
summary(mc)
fviz_mclust(mc, 'classification')
mc$parameters

#results <- results[,2:10]
#write.csv(results, file = outFile)

qprot_result <- read.table('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTout_test', header = TRUE, sep = '\t')
qprot_zeds <- qprot_result[,c(1, 16)]

qprot_dens <- densityMclust(qprot_zeds$Zstatistic)
plot(qprot_dens, what = 'density' )
