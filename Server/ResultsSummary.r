#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(clusterProfiler))
#library(ReactomePA)
#library(data.table)
#library(stringr)
suppressPackageStartupMessages(library(tidyr))
#library(organism, character.only = TRUE)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
suppressPackageStartupMessages(library(organism, character.only = TRUE))

dataSet <- args[1]
#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/PXD004501/Qmodel/PA')

results <- data.frame(matrix(ncol = 9, nrow = 0))
headers <- c('Normalisation',"Threshold", "Go_BP", "Go_MF", "Go_CC",
             "Total", 'DEs', 'DE')
colnames(results) <- headers

norms <- c("None", "Loess-G","RLR-G","VSN-G","TI-G","MedI-G","AI-G","Quantile")

for (norm in norms) {
    setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/PA', sep = ""))
        data <- tryCatch(read.csv(file = paste(dataSet, "_PAout_", norm, ".csv", sep = "")), error = function(e) NULL)
        data <- data[, c(2:9)]
        data$DE <- 'QM'
        results <- rbind(results, data)
}
for (norm in norms) {
    setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Ttest/PA', sep = ""))
    data <- tryCatch(read.csv(file = paste(dataSet, "_T_PAout_", norm, ".csv", sep = "")), error = function(e) NULL)
    data <- data[, c(2:9)]
    data$DE <- 'T'
    results <- rbind(results, data)
}

results <- results[order(-results$Total), ]
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Summary', sep = ""))
write.csv(results, paste(dataSet, '_ResultsSummary.csv', sep = ''))


bestNorm <- results[1,2]
bestThresh <- results[1,3]
bestDE <- results[1,1]

if(bestDE == 'T') {
    #setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/PXD004682/Ttest/Results')
    setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Ttest/Results', sep = ""))
    bestFile <- paste(dataSet, '_', bestNorm, '-normalized_Log2_Tresults.csv', sep = '')
    print(paste('Best result is ', bestFile, ' Threshold: ', bestThresh, ' DE: ', bestDE, sep = ''))
    bestResults <- read.csv(bestFile, header = TRUE)
    DE_prots <- as.vector(subset(bestResults[,2], as.numeric(bestResults$BHpVal) < bestThresh))
    protAcc <- bestResults %>% separate(X0, c("SP", "Uniprot_Acc", "Entrez"), sep = "\\|")
    search_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$BHpVal) < bestThresh))
} else {
    #setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/PXD004501/Qmodel/FDR')
    setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/FDR', sep = ""))
    bestFile <- paste(dataSet, '_', bestNorm, '_FDRresults.csv', sep = '')
    print(paste('Best result is ', bestFile, ' Threshold: ', bestThresh, ' DE: ', bestDE, sep = ''))
    bestResults <- read.csv(bestFile, header = TRUE)
    DE_prots <- as.vector(subset(bestResults$Protein, as.numeric(bestResults$QM_fdr) < bestThresh))
    protAcc <- bestResults %>% separate(Protein, c("SP", "Uniprot_Acc", "Entrez"), sep = "\\|")
    search_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$QM_fdr) < bestThresh))
}
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Summary', sep = ''))
write.csv(DE_prots, file = paste(dataSet, '_DE_prots.csv'))

BG_prots <- protAcc$Uniprot_Acc
DE <- length(search_prots)

if (DE > 0) {
    
    GO_BPresult <- enrichGO(gene = search_prots, 
                            OrgDb = organism, 
                            keyType = 'UNIPROT', 
                            ont = 'BP', 
                            universe = BG_prots, 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'BH')
    if(!is.null(GO_BPresult)) {
        GO_BPsimply <- simplify(GO_BPresult,
                                cutoff = 0.7,
                                by = "p.adjust",
                                select_fun = min)
        GOBP <- data.frame(GO_BPsimply) 
    } 
    GO_MFresult <- enrichGO(gene = search_prots, 
                            OrgDb = organism, 
                            keyType = 'UNIPROT', 
                            ont = 'MF', 
                            universe = BG_prots, 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'BH')
    if(!is.null(GO_MFresult)) {
        GO_MFsimply <- simplify(GO_MFresult,
                                cutoff = 0.7,
                                by = "p.adjust",
                                select_fun = min)
        GOMF <- data.frame(GO_MFsimply)
    } 
    
    GO_CCresult <- enrichGO(gene = search_prots, 
                            OrgDb = organism, 
                            keyType = 'UNIPROT', 
                            ont = 'CC', 
                            universe = BG_prots, 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'BH')
    if(!is.null(GO_CCresult)) {
        GO_CCsimply <- simplify(GO_CCresult,
                                cutoff = 0.7,
                                by = "p.adjust",
                                select_fun = min)
        GOCC <- data.frame(GO_CCsimply)
    }
    totGO <- rbind(GOBP, GOMF, GOCC)
    
} 
write.csv(totGO, file = paste(dataSet, '_PA.csv', sep = ""))
