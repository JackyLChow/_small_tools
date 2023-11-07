signature_calculator <- function(signature_list, # list has names as signature names; signature is character vector
                                 counts = tpm){
  qnt_ <- normalize.quantiles.robust(counts, use.median = T, keep.names = T, n.remove = 0, remove.extreme = "none") # quantile normalize
  qnt_ <- log2(qnt_ + 1) # log2 transform
  
  sig_scores <- list()
  for(i in names(signature_list)){
    genes <- signature_list[[i]]
    missing_genes <- genes[!genes %in% rownames(qnt_)]
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
    if(length(genes) == 0){
    }
    rm(i)
  }
  sig_scores <- data.frame(bind_cols(sig_scores))
  rownames(sig_scores) <- colnames(qnt_)
  return(sig_scores)
}

# # signature scores
# ## calculate scores
# signature_calculator <- function(signature_list, # list has names as signature names; signature is character vector
#                                  expression_table){ # expression table has gene names as rows
#   sig_scores <- list()
#   for(i in names(signature_list)){
#     genes <- signature_list[[i]]
#     missing_genes <- genes[!genes %in% rownames(expression_table)]
#     if(sum(genes %in% rownames(expression_table)) > 1){
#       sig_scores[[i]] <- colMeans(expression_table[genes[genes %in% rownames(expression_table)], ])
#     }
#     if(sum(genes %in% rownames(expression_table)) == 1){
#       sig_scores[[i]] <- expression_table[genes[genes %in% rownames(expression_table)], ]
#     }
#     if(sum(genes %in% rownames(expression_table)) == 0){
#       cat(paste0("genes in ", i, " signature not present in expression table\n"))
#     }
#     rm(i)
#   }
#   sig_scores <- data.frame(bind_cols(sig_scores))
#   rownames(sig_scores) <- colnames(expression_table)
#   return(sig_scores)
# }
