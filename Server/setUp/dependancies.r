options("repos" = c(CRAN = "http://cran.ma.imperial.ac.uk/"))

# install.packages(c("XML",
#                    "unisensR",
#                    "GGIR"))
# 
# library(GGIR)

if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")

BiocManager::install("NormalyzerDE", force = TRUE)
# library(NormalyzerDE)
# BiocManager::install("vsn", force = TRUE)

# BiocManager::install("clusterProfiler")
# organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
#library(organism, character.only = TRUE)

#install.packages('tidyr', repos = 'http://cloud.r-project.org/')
remove.packages('jpeg')
install.packages('jpeg', repos = 'http://cloud.r-project.org/')
library(jpeg)

BiocManager::install("NormalyzerDE", force = TRUE)
library(NormalyzerDE)