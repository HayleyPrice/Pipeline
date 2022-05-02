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

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/QPROTresults')
fileIn <-args[1]
dataSet <- args[2]
fdrThresh <- args[3]
DEmethod <- "QM_FDR"
norm <- "None"
setwd('/mnt/hc-storage/users/hprice/Pipeline/QPROTresults')

fileI <- paste(fileIn, "_FDRresults.csv", sep = '')
#fileOut <- paste (fdrThresh, fileIn, "_PAout.csv", sep = '_')
fileOut <- paste (fdrThresh, fileIn, "_QMfdr_PAout.csv", sep = '_')
# Reads in file and removes top 2 lines
#qprot_out <- read.table(fileIn, header = TRUE)
qprot_out <- read.csv(fileI, header = TRUE)

#qprot_out <- qprot_out[order(qprot_out$fdr), ]
qprot_out <- qprot_out[order(qprot_out$QM_fdr), ]

## create dataframe to store number of significant terms for each threshold
thresholdResults <- data.frame(matrix(ncol = 9, nrow = 0))
headers <- c("DE", "Normalisation", "Threshold", "Go_BP", "Go_MF", "Go_CC",
             "Total", 'DEs')
colnames(thresholdResults) <- headers

#for (fdrThresh in thresholds) {
  
  ## Seperate the protein name so the Unipropt accesion can be used
  protAcc <- qprot_out %>% separate(Protein, c("SP", "Uniprot_Acc", "Entrez"), sep = "\\|")
  #protAcc <- qprot_out %>% separate(Protein, c("num", "Uniprot_Acc"), sep = "\\::")
  
  ## List of background proteins for enrichment analysis
  BG_prots <- protAcc$Uniprot_Acc
  #BG_react <- select(org.Hs.eg.db, BG_prots, "ENTREZID", "UNIPROT")
  
  ## Initialises vector of results
  results <- NULL
  
  # Create list of DE proteins and all proteins
  DE_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$QM_fdr) < fdrThresh))
  #DE_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$fdr) < fdrThresh))           ##QPROT test
  DE <- length(DE_prots)
  NonDE_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$QM_fdr) >= fdrThresh))
  #NonDE_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$fdr) >= fdrThresh))    ##QPROT test
  nonDE <- length(NonDE_prots)
  
  if (DE > 0) {

    GO_BPresult <- enrichGO(gene = DE_prots, 
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
      GOBPnum <- nrow(GO_BPsimply) 
    } else {
      GOBPnum <- 0
    }
    
    GO_MFresult <- enrichGO(gene = DE_prots, 
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
      GOMFnum <- nrow(GO_MFsimply)
    } else {
      GOMFnum <- 0
    }
    
    GO_CCresult <- enrichGO(gene = DE_prots, 
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
      GOCCnum <- nrow(GO_CCsimply)
    }else {
      GOCCnum <- 0
    }
    totGO <- GOBPnum + GOMFnum + GOCCnum

    ## Summary of results
    results <- data.frame(DE = DEmethod, Normalisation = norm, Threshold = fdrThresh, Go_BP = GOBPnum, Go_MF = GOMFnum, 
                          Go_CC = GOCCnum, Total = totGO, DEs = DE)
    
    ## Add to table of results
    thresholdResults <- rbind(thresholdResults, results)
    
  } else {
    print(paste("No DE proteins at threshold: ", fdrThresh, sep = ""))
    flush.console()
    results <- data.frame(DE = DEmethod, Normalisation = norm, Threshold =fdrThresh, Go_BP = 0, Go_MF = 0, 
                          Go_CC = 0, Total = 0, DEs = DE)
    thresholdResults <- rbind(thresholdResults, results)
  }
#}
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/QPROTresults/', dataSet, sep = ""))
write.csv (thresholdResults, file = paste(fileOut, ".csv", sep = ""))
