library(org.Hs.eg.db)
# input is results from GSEA or ORA and the ENTREZID column name
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
