
library(clusterProfiler)
library(ReactomePA)
library(data.table)
library(stringr)
library(tidyr)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/PathwayAnalysis')

## create dataframe to store max number of significant terms for parameter file
paramResults <- data.frame(matrix(ncol = 3, nrow = 0))
#paramHeaders <- c("Parameter", "Max_Terms", "Threshold")
#colnames(paramResults) <- paramHeaders

## Create directories for the results
#fNames <- paste(cutOffs)
#for (name in fNames) {
#  dir.create(name)
#}

fileIn <- "QPROT_FDRtest"                                                               

#Initialise terms
maxTerm <- 0
optThresh <- NULL
paramBest <- NULL

print(paste("Processing file: ", fileIn, sep = ""))
flush.console()
# Reshaping data frames
#######################
# Reads in file and removes top 2 lines
qprot_out <- read.table(fileIn, header = TRUE, sep = "\t")
qprot_out <- qprot_out[,c(1:19, 25)]

# Column names
colnames(qprot_out) <- c("Protein", "N1", "N2", "N3", "N4", "N5", "N6", "N7", 
                         "T1", "T2", "T3", "T4", "T5", "T6", "LogFoldChange", 
                         "Zstatistic", "fdr", "FDRup", "FDRdown", "FDRest")

# Order by pValue
#qprot_out <- qprot_out[order(qprot_out$fdr), ]
qprot_out <- qprot_out[order(qprot_out$FDRest), ]

## Seperate the protein name so the Unipropt accesion can be used
protAcc <- qprot_out %>% separate(Protein, c("SP", "Uniprot_Acc", ".."), sep = "\\|")

## List of background proteins for enrichment analysis
BG_prots <- protAcc$Uniprot_Acc

## create dataframe to store number of significant terms for each threshold
thresholdResults <- data.frame(matrix(ncol = 7, nrow = 0))
headers <- c("Threshold", "Go_BP", "Go_MF", "Go_CC", "React_Path", 
             "Kegg_Path", "Total")
colnames(thresholdResults) <- headers

## Perform enrichment analysis for specified pValue cutoffs
fdrThresh <- 0.005

## Initialises vector of results
results <- NULL

# Create list of DE proteins and all proteins
DE_prots <- as.vector(subset(protAcc$Uniprot_Acc, protAcc$fdr < fdrThresh))
DE <- length(DE_prots)
NonDE_prots <- as.vector(subset(protAcc$Uniprot_Acc, protAcc$fdr >= fdrThresh))
nonDE <- length(NonDE_prots)

if (DE > 0) {
  
  # Get functional annotation chart as R object and save to file
  GOresult <- enrichGO(gene = DE_prots, 
                       OrgDb = 'org.Hs.eg.db', 
                       keyType = 'UNIPROT', 
                       ont = 'ALL', 
                       universe = BG_prots, 
                       pvalueCutoff = 0.05,
                       pAdjustMethod = 'BH')
  
  totGO <- nrow(GOresult)
  
  goBP <- subset(GOresult, GOresult$ONTOLOGY == 'BP')
  goBPNum <- nrow(goBP)
  
  goMF <- subset(GOresult, GOresult$ONTOLOGY =="f")
  goMFNum <- nrow(goMF)
  
  goCC <- subset(FGOresult, GOresultt$ONTOLOGY == "CC")
  goCCNum <- nrow(goCC)
  
  KEGGresult <- enrichKEGG(gene = DE_prots, 
                           organism = "hsa",
                           keyType = 'uniprot',
                           universe = BG_prots,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = 'BH')
  
  keggNum <- nrow(KEGGresult)
  
  REACTresult <- enrichPathway(gene = DE_prots, 
                               organism = "human",
                               universe = BG_prots,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = 'BH')
  
  reactPW <- subset(FuncAnnotChart, Category == "REACTOME_PATHWAY")
  reactPWSigTerms <- subset(reactPW, Benjamini < 0.05)
  reactPWSigNum <- nrow(reactPWSigTerms)
  
  
  ## Summary of results
  results <- data.frame(Threshold = fdrThresh, Go_BP = goBPSigNum, Go_MF = goMFSigNum, 
                        Go_CC = goCCSigNum, React_Path = reactPWSigNum, Kegg_Path = keggPWSigNum, 
                        Total = totSigNum)
  
  #Updates details of best threshold for parameter
  if(totSigNum > maxTerm) {
    maxTerm <- totSigNum
    optThresh <- fdrThresh
  }
  
  ## Add to table of results
  thresholdResults <- rbind(thresholdResults, results)
  
} else {
  results <- data.frame(Threshold = pThresh, Number_Proteins = DEnum, Go_BP = 0, Go_MF = 0, 
                        Go_CC = 0, React_Path = 0, Kegg_Path = 0, 
                        Total = 0)
  thresholdResults <- rbind(thresholdResults, results)
}


#paramBest <- data.frame(Parameter = param, Max_Terms = maxTerm, Threshold = optThresh)
#paramResults <- rbind(paramResults, paramBest)

#write.csv (thresholdResults, file = paste("fdrThreshResults", fileIn, ".csv", sep = ""))
write.table (thresholdResults, file = "PXD004682_QPROT_0_1up.csv", append = TRUE, col.names = FALSE, sep = ",")
#thresholdResults <- read.csv(file = paste("fdrThreshResults", fileIn, ".csv", sep = ""))
#lowThresh <- subset(thresholdResults, thresholdResults$Threshold <= 0.01)
final <- read.table("PXD004682_QPROT_0_1up.csv", header = TRUE, sep = ",")
jpeg(paste("fdrThreshResults", fileIn, ".jpg", sep = ""))
#plot(thresholdResults$Total ~ thresholdResults$Threshold, xlab = "FDR", ylab = "No. Sig Terms")
plot(final$Total ~ final$Threshold, xlab = "FDR", ylab = "No. Sig Terms")
dev.off()


#}

#write.csv (paramResults, file = "fdrParamResults_Burnin.csv")
#thresholdResults <- NULL
