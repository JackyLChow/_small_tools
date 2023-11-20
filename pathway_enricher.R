################################################################################
#
# Pathway enricher
#
################################################################################
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(msigdbr)
library(data.table)
library(dplyr)

# function is to replace ENTREZID with SYMBOL
entrezid_to_symbol <- function(enrichment_results, entrezid_column_name,
                               gene_key = NULL, gene_key_id_column = NULL, gene_key_symbol_column = NULL){
  enrichment_results_ <- enrichment_results
  enrichment_results_[, "symbol"] <- NULL
  for(i in rownames(enrichment_results_)){
    entrezids_ <- unlist(str_split(enrichment_results_[i, entrezid_column_name], "/"))
    if(any(is.null(gene_key, gene_key_id_column, gene_key_symbol_column))){
      symbols_ <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = entrezids_,
                                        columns = "SYMBOL",
                                        keytype = "ENTREZID")[, "SYMBOL"]
    }
    if(!any(is.null(gene_key, gene_key_id_column, gene_key_symbol_column))){
      symbols_ <- gene_key[gene_key[, gene_key_id_column] %in% entrezids_, gene_key_symbol_column]
    }
    enrichment_results_[i, "symbol"] <- paste(symbols_, collapse = "/")
  }
  return(enrichment_results_)
}

# input is results from differential gene expression, requires ENTREZID column and to identify ranking column
pathway_enricher <- function(dge_results,
                             ranking_column = "logFC", entrez_id_column = "ENTREZID",
                             pval_cut = 0.1, convert_to_symbol = T){
  if(!is.numeric(dge_results[, entrez_id_column])){
    cat("\nGene ID column must be Entrez ID and numeric\n")
  }
  
  # load hallmark pathways
  msig_h_ <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::rename(ont = gs_name, gene = entrez_gene)
  
  # object to be returned
  pe_output <- list() 
  
  # make ranked gene list
  ## extract ordered gene list by ranking_column
  ordered_genes_ <- dge_results
  ordered_genes_ <- data.frame(ordered_genes_[, entrez_id_column], ordered_genes_[, ranking_column])
  names(ordered_genes_) <- c("ENTREZID", "ranking_column")
  ordered_genes_ <- ordered_genes_[order(ordered_genes_$ranking_column, decreasing = T), ]

  ## make ranked list
  ranked_list_ <- ordered_genes_[, "ranking_column"]
  names(ranked_list_) <- as.character(ordered_genes_[, "ENTREZID"])
  
  # run pathway enrichments
  ## reactome pathway enrichment
  set.seed(415); pe_output[["Reactome_GSEA"]] <- gsePathway(ranked_list_,
                                                            maxGSSize = 500,
                                                            pvalueCutoff = pval_cut) %>% data.frame()
  
  ## KEGG pathway enrichment
  set.seed(415); pe_output[["KEGG_GSEA"]] <- gseKEGG(ranked_list_,
                                                     organism = "hsa",
                                                     pvalueCutoff = pval_cut) %>% data.frame()
  
  ## GO pathway enrichment
  set.seed(415); pe_output[["GO_GSEA"]] <- gseGO(ranked_list_,
                                                 OrgDb = org.Hs.eg.db,
                                                 pvalueCutoff = pval_cut) %>% data.frame()
  
  ## Hallmark pathway enrichment
  set.seed(415); pe_output[["Hallmark_GSEA"]] <- GSEA(ranked_list_, TERM2GENE = msig_h_, scoreType = "pos") %>% data.frame()
  
  # polish
  pe_output <- data.frame(rbindlist(pe_output[sapply(pe_output, nrow) > 0], idcol = "geneset_database"))
  if(convert_to_symbol == T){
    pe_output <- entrezid_to_symbol(pe_output, "core_enrichment")
  }
  
  return(pe_output)
  rm(ordered_genes_, ranked_list_, pe_output, key_, msig_h_)
}
