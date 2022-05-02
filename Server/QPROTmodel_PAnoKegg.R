#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(clusterProfiler))
#library(ReactomePA)
#library(data.table)
#library(stringr)
suppressPackageStartupMessages(library(tidyr))
#library(organism, character.only = TRUE)

# SET THE DESIRED ORGANISM HERE
#organism = "org.Mm.eg.db"
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
suppressPackageStartupMessages(library(organism, character.only = TRUE))

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/PXD004501/Qmodel/FDR')
norm <-args[1]
dataSet <- args[2]
DEmethod <- "QPROTmodel"
#norm <- "VSN-G"
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/FDR', sep = ""))

fileIn <- paste(dataSet, "_", norm, "_FDRresults.csv", sep = '')
fileOut <- paste (dataSet, "_PAout_", norm, sep = '')
# Reads in file and removes top 2 lines
qprot_out <- read.csv(fileIn, header = FALSE)
cols <- qprot_out[1,]
qprot_out <- qprot_out[-1,]

index1 <- grep("1", cols)
ind1 <- length(index1)
index2 <- grep("2", cols)
ind2 <- length(index2)

if(ind1 < ind2){
  index <- ind2
} else {
  index <- ind1
}

s <- LETTERS[seq(1,index)]

for (i1 in 1:ind1) {
  j1 <- index1[i1]
  cols[1,j1] <- paste(cols[1,j1], s[i1], sep = '')
}
for (i2 in 1:ind2) {
  j2 <- index2[i2]
  cols[1,j2] <- paste(cols[1,j2], s[i2], sep = '')
}

# Column names
colnames(qprot_out) <- cols
qprot_out <- qprot_out[order(qprot_out$QM_fdr), ]

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
  #protAcc <- qprot_out %>% separate(Protein, "Uniprot_Acc", sep = ";")
  #protAcc <- qprot_out %>% separate(Protein, c("num", "Uniprot_Acc"), sep = "\\::")
  
  ## List of background proteins for enrichment analysis
  BG_prots <- protAcc$Uniprot_Acc
  #BG_react <- select(org.Hs.eg.db, BG_prots, "ENTREZID", "UNIPROT")
  
  ## Initialises vector of results
  results <- NULL
  #QM_fdr <-0.000001
  
  # Create list of DE proteins and all proteins
  DE_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$QM_fdr) < fdrThresh))
  DE <- length(DE_prots)
  NonDE_prots <- as.vector(subset(protAcc$Uniprot_Acc, as.numeric(protAcc$QM_fdr) >= fdrThresh))
  nonDE <- length(NonDE_prots)
  
  if (DE > 0) {
    
    print(paste(norm, ' threshold ', fdrThresh, "- Number of DEs: ", DE, sep = ""))
    flush.console()
    
    #print("Performing GO analysis....")
    flush.console()
    # Get functional annotation chart as R object and save to file
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
    print(paste(norm, ' threshold ', fdrThresh, " - Number of GO terms: ",totGO, sep = ""))
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
}
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/PA', sep = ""))
write.csv (thresholdResults, file = paste(fileOut, ".csv", sep = ""))
