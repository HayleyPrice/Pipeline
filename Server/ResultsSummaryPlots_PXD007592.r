#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(clusterProfiler))
library(ggplot2)
#library(ReactomePA)
#library(data.table)
#library(stringr)
suppressPackageStartupMessages(library(tidyr))
library(ggh4x)
library(wesanderson)
#library(organism, character.only = TRUE)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
suppressPackageStartupMessages(library(organism, character.only = TRUE))

#dataSet <- args[1]
#dataSet <- "PXD004501"
#dataSet <- "PXD004682"
dataSet <- "PXD007592"
    
results <- data.frame(matrix(ncol = 9, nrow = 0))
headers <- c('Normalisation',"Threshold", "Go_BP", "Go_MF", "Go_CC",
             "Total", 'DEs', 'DE')
colnames(results) <- headers

norms <- c("None", "AI-G", "MedI-G", "TI-G", "RLR-G", "Loess-G","VSN-G","Quantile")
labels <- c("Log2 transformed", 
            "Average intensity normalisation",
            "Median intensity normalisation", 
            "Total intensity normalisation",
            "Robust linear regression normalisation",
            "Loess normalisation",
            "Variance stabilisation normalisation",
            "Quantile normalisation")

for (norm in norms) {
    #setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/PA', sep = ""))
    setwd(paste('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/', dataSet, '/Qmodel/PA', sep = ""))
        data <- tryCatch(read.csv(file = paste(dataSet, "_PAout_", norm, ".csv", sep = "")), error = function(e) NULL)
        data <- data[, c(2:9)]
        data$DE <- 'BT'
        results <- rbind(results, data)
}
for (norm in norms) {
    #setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Ttest/PA', sep = ""))
    setwd(paste('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/', dataSet, '/Ttest/PA', sep = ""))
    data <- tryCatch(read.csv(file = paste(dataSet, "_T_PAout_", norm, ".csv", sep = "")), error = function(e) NULL)
    data <- data[, c(2:9)]
    data$DE <- 'T'
    results <- rbind(results, data)
}

results <- results[order(-results$Total), ]
#results <- results[order(-results$DEs), ]
#results <- results[-c(1,2), ]

#setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Summary', sep = ""))
#write.csv(results, paste(dataSet, '_ResultsSummary.csv', sep = ''))
resultsP <- results[,c(1:3,7,8)]
resultsP$DataSet <- dataSet
r1.melt <- melt(resultsP[,c("Threshold", "DEs", "Total", "Normalisation", "DataSet", "DE")], id.vars = c(1,2,3, 4, 5, 6))
r1.melt$Normalisation <- factor(r1.melt$Normalisation, levels = norms, labels = labels)
r1.melt$DE <- factor(r1.melt$DE, levels = c("BT", "T"), labels = c("BayT", "t-Test"))

p1 <- ggplot(data = r1.melt, aes(x = Threshold)) +
    geom_point(aes(y = Total, color = DE)) +
    geom_line(aes(y = DEs/10, color = DE)) +
    facet_grid(. ~ Normalisation, labeller = label_wrap_gen(width = 19), scales = "free") +
    #scale_x_log10(name = 'Threshold (log scale)', breaks = c( 0.001, 0.01, 0.1), labels = c( 0.001, 0.01, 0.1)) + 
    scale_x_log10(name = 'Threshold (log scale)') +
    scale_y_continuous(name = 'Number of significant terms',
                       sec.axis = sec_axis(~.*10, name="Number of DE proteins")) +
    theme_bw() + 
    theme(strip.text.y = element_text(size = 16, face = 'bold'),
          strip.text.x = element_text(size = 12, face = 'bold'),
          axis.title.x = element_text(size = 16, face = 'bold'),
          axis.text.x = element_text(size = 11, face = 'bold'),
          axis.title.y = element_text(size = 16, face = 'bold'),
          axis.text.y = element_text(size = 11, face = 'bold'),
          legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14))

p1 <- p1 +  scale_colour_manual(values=c("darkorange3", "cadetblue"), 
                                name="Differential Expression Method",
                                labels = c("BayT", "t-Test"))
setwd(paste('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/', dataSet, '/Summary', sep = ""))

png("TermsPlot.png", height = 750, width = 1400)
print(p1)
dev.off()

bestNorm <- results[1,2]
bestThresh <- results[1,3]
bestDE <- results[1,1]

if(bestDE == 'T') {
    setwd(paste('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/', dataSet, '/Ttest/Results', sep = ""))
    #setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Ttest/Results', sep = ""))
    bestFile <- paste(dataSet, '_', bestNorm, '-normalized_Log2_Tresults.csv', sep = '')
    print(paste('Best result is ', bestFile, ' Threshold: ', bestThresh, ' DE: ', bestDE, sep = ''))
    bestResults <- read.csv(bestFile, header = TRUE)
    DE_prots <- as.vector(subset(bestResults[,2], as.numeric(bestResults$BHpVal) < bestThresh))
    protAcc <- bestResults %>% separate(X0, c("SP", "Uniprot_Acc", "Entrez"), sep = "\\|")
    search_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$BHpVal) < bestThresh))
} else {
    setwd(paste('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/', dataSet, '/Qmodel/FDR', sep = ""))
    #setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/FDR', sep = ""))
    bestFile <- paste(dataSet, '_', bestNorm, '_FDRresults.csv', sep = '')
    print(paste('Best result is ', bestFile, ' Threshold: ', bestThresh, ' DE: ', bestDE, sep = ''))
    bestResults <- read.csv(bestFile, header = TRUE)
    DE_prots <- as.vector(subset(bestResults$Protein, as.numeric(bestResults$QM_fdr) < bestThresh))
    protAcc <- bestResults %>% separate(Protein, c("SP", "Uniprot_Acc", "Entrez"), sep = "\\|")
    search_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$QM_fdr) < bestThresh))
}
#setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Summary', sep = ''))
setwd(paste('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/', dataSet, '/Summary', sep = ""))
#write.csv(DE_prots, file = paste(dataSet, '_DE_prots.csv'))

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
        png("GO.png")
        print(barplot(GO_BPsimply))
        dev.off()
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
    png("MF.png")
    print(barplot(GO_MFsimply))
    dev.off()
    
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
    png("CC.png")
    print(barplot(GO_CCsimply))
    dev.off()
    
    totGO <- rbind(GOBP, GOMF, GOCC)
    
} 
write.csv(totGO, file = paste(dataSet, '_PA.csv', sep = ""))
