# Explore gene expression across TCGA projects

### load data ------------------------------------------------------------------
counts <- readRDS("~/Documents/BFX_proj/_data_public/TCGA_PanCancer/_processed/tcga_pancancer_rsem_log_tf_FZ01.rds")
metadata <- readRDS("~/Documents/BFX_proj/_data_public/TCGA_PanCancer/_processed/tcga_pancancer_metadata_FZ01.rds")
rownames(metadata) <- metadata$SAMPLE_ID

tcga_expression <- function(counts, gene){
  df_ <- data.frame(metadata, expression = counts[gene, ])
  cancer_type_acronym <- c()
  mean_expression <- c()
  for(i in unique(metadata$CANCER_TYPE_ACRONYM)){
    cancer_type_acronym <- c(cancer_type_acronym, i)
    mean_expression <- c(mean_expression, mean(df_[df_$CANCER_TYPE_ACRONYM == i, "expression"]))
  }
  df_$CANCER_TYPE_ACRONYM <- factor(df_$CANCER_TYPE_ACRONYM,
                                    levels = cancer_type_acronym[order(mean_expression, decreasing = T)])
  ggplot(df_, aes(CANCER_TYPE_ACRONYM, expression)) +
    ylab(gene) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.y = element_text(face = "italic"),
          axis.title.x = element_blank())
}

tcga_expression(counts = counts, gene = "IFNG")
