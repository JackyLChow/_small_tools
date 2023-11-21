library(preprocessCore)

signature_calculator <- function(signature_list, # list has names as signature names; signature is character vector of genes
                                 counts = log2_tpm_plus_one){
  qnt_ <- normalize.quantiles.robust(counts, use.median = T, keep.names = T, n.remove = 0, remove.extreme = "none") # quantile normalize
  qnt_ <- log2(qnt_ + 1) # log2 transform
  
  sig_scores <- list()
  for(i in names(signature_list)){
    genes <- signature_list[[i]]
    missing_genes <- genes[!genes %in% rownames(qnt_)]
    if(length(missing_genes == length(genes)){
      cat(paste("all genes missing for", i, "\n"))
    }
    if(length(missing_genes) > 0){
      cat(paste("missing genes for ", i, ":", missing_genes, "\n"))
    }
    if(length(missing_genes) == 0){
      cat(paste("all genes present for", i, "\n"))
    }
    genes <- genes[genes %in% rownames(qnt_)]
    if(length(genes) > 1){
      sig_scores[[i]] <- colMeans(qnt_[genes, ])
    }
    if(length(genes) == 1){
      sig_scores[[i]] <- qnt_[genes, ]
    }
    rm(i)
  }
  sig_scores <- data.frame(bind_cols(sig_scores))
  rownames(sig_scores) <- colnames(qnt_)
  return(sig_scores)
}