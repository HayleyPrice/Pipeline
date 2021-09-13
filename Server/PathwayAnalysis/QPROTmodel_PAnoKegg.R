args = commandArgs(trailingOnly=TRUE)

library(clusterProfiler)
#library(ReactomePA)
#library(data.table)
#library(stringr)
library(tidyr)
#library(organism, character.only = TRUE)
library(tidyr)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/PathwayAnalysis')
norm <-args[2]
DE <- "QPROTmodel"
#norm <- "RLR-G"

fileIn <- paste(norm, "_DEresults.csv", sep = '')
#fileOut <- paste ("PAout_", norm, sep = '')
fileOut <- paste ("results/PAout_", norm, sep = '')
# Reads in file and removes top 2 lines
qprot_out <- read.csv(fileIn, header = TRUE)
qprot_out <- qprot_out[,c(3:20)]

# Column names
colnames(qprot_out) <- c("Protein", "N1", "N2", "N3", "N4", "N5", "N6", "N7", 
                         "T1", "T2", "T3", "T4", "T5", "T6", "LogFoldChange", 
                         "Zstatistic", "pVal", "BHpVal")

thresholds <- c(0.000001, 0.000002, 0.000003, 0.000004, 0.000005, 0.000006, 0.000007, 0.000008, 0.000009,
                0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008, 0.00009,
                0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009,
                0.001,  0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009,
                0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)

## create dataframe to store number of significant terms for each threshold
thresholdResults <- data.frame(matrix(ncol = 9, nrow = 0))
headers <- c("DE", "Normalisation", "Threshold", "Go_BP", "Go_MF", "Go_CC",
             "Total", 'DEs')
colnames(thresholdResults) <- headers

for (fdrThresh in thresholds) {

  ## Seperate the protein name so the Unipropt accesion can be used
  protAcc <- qprot_out %>% separate(Protein, c("SP", "Uniprot_Acc", "Entrez"), sep = "\\|")
  
  ## List of background proteins for enrichment analysis
  BG_prots <- protAcc$Uniprot_Acc
  BG_react <- select(org.Hs.eg.db, BG_prots, "ENTREZID", "UNIPROT")

  ## Initialises vector of results
  results <- NULL
  
  # Create list of DE proteins and all proteins
  DE_prots <- as.vector(subset(protAcc$Uniprot_Acc, protAcc$BHpVal < fdrThresh))
  DE <- length(DE_prots)
  NonDE_prots <- as.vector(subset(protAcc$Uniprot_Acc, protAcc$BHpVal >= fdrThresh))
  nonDE <- length(NonDE_prots)
  
  if (DE > 0) {
    
    print(paste("Number of DEs: ", DE, sep = ""))
    flush.console()
    
    print("Performing GO analysis....")
    flush.console()
    # Get functional annotation chart as R object and save to file
    GO_BPresult <- enrichGO(gene = DE_prots, 
                            OrgDb = 'org.Hs.eg.db', 
                            keyType = 'UNIPROT', 
                            ont = 'BP', 
                            universe = BG_prots, 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'BH')
    
    GO_BPsimply <- simplify(GO_BPresult,
                            cutoff = 0.7,
                            by = "p.adjust",
                            select_fun = min)
    GOBPnum <- nrow(GO_BPsimply) 
    
    GO_MFresult <- enrichGO(gene = DE_prots, 
                            OrgDb = 'org.Hs.eg.db', 
                            keyType = 'UNIPROT', 
                            ont = 'MF', 
                            universe = BG_prots, 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'BH')
    GO_MFsimply <- simplify(GO_MFresult,
                            cutoff = 0.7,
                            by = "p.adjust",
                            select_fun = min)
    GOMFnum <- nrow(GO_MFsimply)
    
    GO_CCresult <- enrichGO(gene = DE_prots, 
                            OrgDb = 'org.Hs.eg.db', 
                            keyType = 'UNIPROT', 
                            ont = 'CC', 
                            universe = BG_prots, 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'BH')
    GO_CCsimply <- simplify(GO_CCresult,
                            cutoff = 0.7,
                            by = "p.adjust",
                            select_fun = min)
    GOCCnum <- nrow(GO_CCsimply)
    
    totGO <- GOBPnum + GOMFnum + GOCCnum
    print(paste("Number of GO terms: ", totGO, sep = ""))
    flush.console()
    
    # print("Performing KEGG analysis")
    # flush.console()
    # KEGGresult <- enrichKEGG(gene = DE_prots, 
    #                          organism = "hsa",
    #                          keyType = 'uniprot',
    #                          universe = BG_prots,
    #                          pvalueCutoff = 0.05,
    #                          pAdjustMethod = 'BH')
    # KEGGnum <- nrow(KEGGresult)
    # print(paste("Number of KEGG terms: ", KEGGnum , sep = ""))
    # flush.console()
    # 
    # totNum <- totGO + KEGGnum
    
    #   ## Summary of results
    #   results <- data.frame(Threshold = fdrThresh, Go_BP = GOBPnum, Go_MF = GOMFnum, 
    #                         Go_CC = GOCCnum,  Kegg_Path = KEGGnum , Total = totNum, DEs = DE)
    #   
    #   ## Add to table of results
    #   thresholdResults <- rbind(thresholdResults, results)
    #   
    # } else {
    #   results <- data.frame(Threshold =fdrThresh, Go_BP = 0, Go_MF = 0, 
    #                         Go_CC = 0, Kegg_Path = 0, Total = 0, DEs = DE)
    #   thresholdResults <- rbind(thresholdResults, results)
    ## Summary of results
    results <- data.frame(DE = DE, Normalisation = norm, Threshold = fdrThresh, Go_BP = GOBPnum, Go_MF = GOMFnum, 
                          Go_CC = GOCCnum, Total = totGO, DEs = DE)
    
    ## Add to table of results
    thresholdResults <- rbind(thresholdResults, results)
    
  } else {
    print(paste("No DE proteins at threshold: ", fdrThresh, sep = ""))
    flush.console()
    results <- data.frame(DE = DE, Normalisation = norm, Threshold =fdrThresh, Go_BP = 0, Go_MF = 0, 
                          Go_CC = 0, Total = 0, DEs = DE)
    thresholdResults <- rbind(thresholdResults, results)
  }
}
write.csv (thresholdResults, file = paste(fileOut, ".csv", sep = ""))
