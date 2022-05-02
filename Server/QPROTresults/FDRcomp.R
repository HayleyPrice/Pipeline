library(ggplot2)

setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/QPROTresults')

datasets <- c('PXD004501', 'PXD004682', 'PXD007592')
results <- data.frame()
cols <- c('fdr', 'DE', 'Data')

for(data in datasets) {
    
    QPROT <- read.table(paste(data, '_1up_Normalyser_qprot_qprot_fdr', sep = ''), 
                        header=TRUE, sep="\t")
    QPROTfdr <- as.data.frame(sort(QPROT$fdr))
    colnames(QPROTfdr) <- 'fdr'
    QPROTfdr$Obs <- 1:nrow(QPROTfdr)
    QPROTfdr$DE <- 'QPROT'
    QPROTfdr$Data <- data
    
    QM <- read.table(paste(data, '_1up_Normalyser_qprot_qprot_fdr','_FDRresults.csv', sep = ''), 
                     header=TRUE, sep=",")
    QMfdr <- as.data.frame(sort(QM$QM_fdr))
    colnames(QMfdr) <- 'fdr'
    QMfdr$Obs <- 1:nrow(QMfdr)
    QMfdr$DE <- 'QM'
    QMfdr$Data <- data
    
    results <- rbind(results, QPROTfdr, QMfdr)
}

results$DE <- factor(results$DE, levels = c('QPROT', 'QM'), labels = c('QPROT', 'QPROT\nModel'))
results$Data <- factor(results$Data, levels = datasets)

p <- ggplot(results, aes(x = DE, y = fdr)) +
    geom_violin() +
    facet_grid(cols = vars(Data)) +
    scale_x_discrete(name = '') +
    scale_y_continuous(name = 'False discovery rate') +
    theme_bw() +
    theme(strip.text.x = element_text(size = 14, face = 'bold'),
          axis.text.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 14, face = 'bold'),
          axis.text.y = element_text(size = 12, face = 'bold'))

png('VioloinPlot.png', height = 400, width = 700)
print(p)
dev.off()

q <- ggplot(results, aes(x = DE, y = fdr)) +
    geom_boxplot() +
    facet_grid(cols = vars(Data)) +
    scale_x_discrete(name = '') +
    scale_y_continuous(name = 'False discovery rate') +
    theme_bw() +
    theme(strip.text.x = element_text(size = 14, face = 'bold'),
          axis.text.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 14, face = 'bold'),
          axis.text.y = element_text(size = 12, face = 'bold'))

png('BoxPlot.png', height = 500, width = 700)
print(q)
dev.off()

r <- ggplot(results, aes(x = Obs, y = fdr)) +
    geom_point(aes(color = DE)) +
    facet_grid(cols = vars(Data)) +
    scale_y_continuous(name = 'FDR') +
    scale_x_continuous(name = 'Observations') +
    theme_bw() +
    labs(color='FDR method') +
    theme(strip.text.x = element_text(size = 14, face = 'bold'),
          axis.text.x = element_text(size = 12, face = 'bold'),
          axis.title.x = element_text(size = 14, face = 'bold'),
          axis.title.y = element_text(size = 14, face = 'bold'),
          axis.text.y = element_text(size = 12, face = 'bold'),
          legend.position = 'bottom',
          legend.title = element_text(size = 16, face = 'bold'),
          legend.text = element_text(size = 14, face = 'bold'))

png('FDRComp.png', height = 400, width = 700)
print(r)
dev.off()

