#Imports

libs <- c("tidyverse", "ggVennDiagram", "BiocManager",
          "DESeq2","ggplot2","preprocessCore","fgsea")
for (package in libs) {
  suppressPackageStartupMessages(require(package, 
                                         quietly = T, 
                                         character.only = T))
  require(package, character.only = T)
}

#Loading the data, cleaning the data (removing X.1 and X.2 columns) and filtering the zero counts genes along with dropping of the duplicates

counts <- read.csv('updated_counts.csv')
colData <- read.csv('sample_info.csv')
counts = counts[-c(28,29)]
filtered_counts <- counts[rowSums(counts[-c(1, 2,26,27)]) > 0, ]
filtered_counts = filtered_counts[!duplicated(filtered_counts[c(-1)]),]

#Assigning the rowData and the colData

samples <- colData$Sample

rowData <- filtered_counts["SPU.name"]
uniprothit <- filtered_counts["uniprothit"]

colnames(filtered_counts) <- gsub("^sample.", "", colnames(filtered_counts))

colData$Sample <- gsub("-",".",samples)

#Function to perform the pairwise DE analysis using DESeq2 ande the quantile normalization

run_de_pairwise <- function(g1, g2){
  subset_meta_data <- colData[colData$Group==g1 | colData$Group==g2,]
  subset_count_data <- filtered_counts[-c(1,2,26,27)][subset_meta_data$Sample]
  subset_count_data_normalized <- normalize.quantiles(as.matrix(subset_count_data))
  subset_count_data_normalized <- as.matrix(apply(subset_count_data_normalized, 2, function(x) as.integer(x)))
  colnames(subset_count_data_normalized) <- colnames(subset_count_data)
  
  se <- SummarizedExperiment(assays=list(subset_count_data_normalized),colData=subset_meta_data,rowData=rowData)
  se$Group <- factor(se$Group)
  dds <- DESeqDataSet(se, design = ~Group)
  sizeFactors(dds) <- rep(1, ncol(dds))
  dds <- DESeq(dds)
  
  deseq_results <- data.frame(results(dds))
  deseq_results["Genes"] = rowData
  deseq_results["uniprothit"] = uniprothit
  
  
  padj_threshold <- 0.05
  
  deseq_results <-deseq_results[complete.cases(deseq_results), ]
  
  deseq_results$volc_plot_status <- NA 

  deseq_results$volc_plot_status[deseq_results$padj < padj_threshold & deseq_results$log2FoldChange > 0] <- 'UP'
  
  deseq_results$volc_plot_status[deseq_results$padj < padj_threshold & deseq_results$log2FoldChange < 0] <- 'DOWN'
  
  

  deseq_results$volc_plot_status[deseq_results$padj > padj_threshold] <- 'NS'
  volcano_plot <- ggplot(deseq_results, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = volc_plot_status), size = 2, alpha = 0.6) +
    scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "gray")) +
    labs(x = "Log2 Fold Change",
         y = "-log10(padj)",
         title =  paste("Volcano Plot -", g1, "vs", g2)) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white")) 
  
  ggsave(paste("Volcano Plot -", g1, "vs", g2, ".png"), plot = volcano_plot, 
         width = 8, height = 6, units = "in", bg = "white")
  return(deseq_results)
}

