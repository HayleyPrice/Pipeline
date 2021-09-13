library(mclust)
library(factoextra)
library(MASS)
library(ggplot2)
library(FDRestimation)
library(fdrtool)
library(qvalue)

#dataSet <- args[1]
#comp <- args[2]

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/DE/Bayes/QPROTmodel/LocalScripts/test') 

outFile <- "PXD004682_QPROTmodel_testResult.csv"

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

#qprot_qval <- qvalue(qprot_result$pees)
set.seed(1977)
n <- nrow(qprot_result)
P <- 100
variable <- qprot_result$Zstatistic

permSamples <- matrix(0, nrow = n, ncol = P)

for(i in 1:P){
    permSamples[,i] <- sample(variable, size= n, replace=FALSE)
}

qprot_result$fdr_pval <- empPvals(stat = qprot_result$Zstatistic, stat0 = permSamples, pool = FALSE)
pi0 <- get.pi0(pvalues = qprot_result$fdr_pval, set.pi0 = 1, zvalues = qprot_result$Zstatistic, estim.method = "storey")
qprot_est <- p.fdr(pvalues = qprot_result$fdr_pval, zvalues = qprot_result$Zstatistic, set.pi0 = pi0)
qprot_result$fdr_est <- qprot_est$fdrs

qprot_result$pees <- 2*pnorm(-abs(qprot_result$Zstatistic))
pi0_2 <- get.pi0(pvalues = qprot_result$pees, set.pi0 = 1, zvalues = qprot_result$Zstatistic, estim.method = "storey")
qprot_est_2 <- p.fdr(pvalues = qprot_result$pees, zvalues = qprot_result$Zstatistic, set.pi0 = pi0_2)
qprot_result$fdr_est_2 <- qprot_est_2$fdrs


p <- ggplot(qprot_result) +
    geom_density(aes(x = fdr), colour = 'red') +
    geom_density(aes(x = fdr_est_2), colour = 'pink') +
    geom_density(aes(x = fdr_est), colour = 'blue') 
    
p

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/PathwayAnalysis')
write.table(qprot_result, "QPROT_FDRtest", sep = "\t")



