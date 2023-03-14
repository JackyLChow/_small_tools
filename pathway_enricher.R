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
entrezid_to_symbol <- function(enrichment_results, entrezid_column_name){
  enrichment_results_ <- enrichment_results
  enrichment_results_[, "symbol"] <- NULL
  for(i in rownames(enrichment_results_)){
    entrezids_ <- unlist(str_split(enrichment_results_[i, entrezid_column_name], "/"))
    symbols_ <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = entrezids_,
                                      columns = "SYMBOL",
                                      keytype = "ENTREZID")[, "SYMBOL"]
    enrichment_results_[i, "symbol"] <- paste(symbols_, collapse = "/")
  }
  return(enrichment_results_)
}

# input is results from differential gene expression, requires gene symbol as rownames and to identify ranking column

pathway_enricher <- function(dge_results, ranking_column = "logFC", pval_cut = 0.1, convert_to_symbol = T){
  # object to be returned
  pe_output <- list() 
  
  # make ranked gene list
  ## extract ordered gene list by ranking_column
  ordered_genes_ <- dge_results
  ordered_genes_ <- data.frame(Symbol = rownames(ordered_genes_), ranking_column = ordered_genes_[, ranking_column])
  ordered_genes_ <- ordered_genes_[order(ordered_genes_$ranking_column, decreasing = T), ]
  
  ## generate ENTREZID key from synonym table
  key_ <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = ordered_genes_$Symbol,
                                columns = c("ENTREZID"),
                                keytype = "SYMBOL")
  key_ <- key_[!duplicated(key_$SYMBOL), ] # remove duplicated symbols
  
  ## load hallmark pathways
  msig_h_ <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::rename(ont = gs_name, gene = entrez_gene)
  
  ## assign ENTREZID to gene
  ordered_genes_ <- left_join(ordered_genes_, key_, by = c("Symbol" = "SYMBOL"))
  
  ## filter genes with duplicated and missing ENTREZID
  ordered_genes_ <- ordered_genes_[!is.na(ordered_genes_$ENTREZID), ]
  ordered_genes_ <- ordered_genes_[!duplicated(ordered_genes_$ENTREZID), ]
  
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
