pathway_enricher <- function(dge_results, ranking_column = "logFC"){
  # object to be returned
  pe_out <- list() 
  
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
  set.seed(415); reactomeGSEA_ <- gsePathway(ranked_list_,
                                             maxGSSize = 500,
                                             pvalueCutoff = 0.1) %>% data.frame()
  pe_out[["Reactome_GSEA"]] <- entrezid_to_symbol(reactomeGSEA_)
  
  ## KEGG pathway enrichment
  set.seed(415); keggGSEA_ <- gseKEGG(ranked_list_,
                                      organism = "hsa",
                                      pvalueCutoff = 0.1) %>% data.frame()
  pe_out[["KEGG_GSEA"]] <- entrezid_to_symbol(keggGSEA_)
  
  ## GO pathway enrichment
  set.seed(415); goGSEA_ <- gseGO(ranked_list_,
                                  OrgDb = org.Hs.eg.db,
                                  pvalueCutoff = 0.1) %>% data.frame()
  pe_out[["GO_GSEA"]] <- entrezid_to_symbol(goGSEA_)
  
  ## Hallmark pathway enrichment
  set.seed(415); hallmarkGSEA_ <- GSEA(ranked_list_, TERM2GENE = msig_h_, scoreType = "pos") %>% data.frame()
  pe_out[["Hallmark_GSEA"]] <- entrezid_to_symbol(hallmarkGSEA_)
  
  pe_out <- rbindlist(pe_out, idcol = "GSEA")
  
  return(pe_out)
}
