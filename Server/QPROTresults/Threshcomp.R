library(ggplot2)
library(scales)

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/QPROTresults')

datasets <- c('PXD004501', 'PXD004682', 'PXD007592')
results <- data.frame()
#cols <- c('fdr', 'DE', 'Data')

for(data in datasets) {
    
    QPROT <- read.csv(paste(data, '_DE.QP_FDR.QM/ResultsSummary_QPROT_test.csv', sep = ''), header=TRUE)
    QPROT <- QPROT[,c(-1,-2)]
    QPROT$Data <- data
    
    QM <- read.table(paste('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/', 
                           data, '/Qmodel/PA/', data, '_PAout_None.csv', sep = ''), 
                     header=TRUE, sep=",")
    QM <- QM[,-1]
    QM$DE <- 'QM'
    QM$Data <- data

    results <- rbind(results, QPROT, QM)
}

results$DE <- factor(results$DE, levels = c('QM_FDR', 'QM'), labels = c('QPROT', 'QPROT\nModel'))
results$Data <- factor(results$Data, levels = datasets)
results <- results[order(-results$Total),]
PXD004501 <- results[results$Data == 'PXD004501',]
PXD004682 <- results[results$Data == 'PXD004682',]

r <- ggplot(results, aes(x = Threshold, y = Total)) +
    geom_point(aes(color = DE)) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(name = 'Total number of terms') +
    facet_grid(rows = vars(Data)) +
    theme_bw() +
    labs(color='Differential expression method') +
    theme(strip.text.y = element_text(size = 14, face = 'bold'),
          axis.text.x = element_text(size = 12, face = 'bold'),
          axis.title.x = element_text(size = 14, face = 'bold'),
          axis.title.y = element_text(size = 14, face = 'bold'),
          axis.text.y = element_text(size = 12, face = 'bold'),
          legend.position = 'bottom',
          legend.title = element_text(size = 16, face = 'bold'),
          legend.text = element_text(size = 14, face = 'bold'))

png('ThreshComp.png', height = 500, width = 700)
print(r)
dev.off()

