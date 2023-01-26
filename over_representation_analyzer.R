################################################################################
#
# Over representation analyzer
#
################################################################################

library(org.Hs.eg.db)

# input is results from DGE in RNAseq pipeline
data <- readRDS("~/Documents/BFX_proj/RNAseq_pipeline/_output/Chow_PNAS_2020/differential_gene_expression_results.RDS")
data <- data$differential_gene_expression_resutls

# separate probe genes
probe <- rownames(data)[data$logFC > 1 & data$adj.P.Val < 0.05]
universe <- rownames(data)

# function to convert from gene symbols to ENTREZID
entrez_ider <- function(genes){
  return(AnnotationDbi::select(org.Hs.eg.db,
                               keys = genes,
                               columns = "ENTREZID",
                               keytype = "SYMBOL")[, "ENTREZID"])
}

# target for output
over_representation_analyzer <- function(probe, universe){
  output <- list()
  probe_ <- entrez_ider(probe)
  universe_ <- entrez_ider(universe)
  output[["Reactome_OR"]] <- data.frame(enrichPathway(gene         = probe_,
                                                      universe     = universe_,
                                                      pvalueCutoff = 0.1))
  output[["KEGG_OR"]] <- data.frame(enrichKEGG(gene         = probe_,
                                               universe     = universe_,
                                               organism     = 'hsa',
                                               pvalueCutoff = 0.1))
  output[["GO_OR"]] <- data.frame(enrichGO(gene          = probe_,
                                           universe      = universe_,
                                           OrgDb         = org.Hs.eg.db,
                                           pvalueCutoff  = 0.1))
  return(output)
}

#foo <- over_representation_analyzer(probe, universe)
