# signature scores
## calculate scores
signature_calculator <- function(signature_list, expression_table){
  sig_scores <- list()
  for(i in names(signature_list)){
    genes <- signature_list[[i]]
    if(sum(genes %in% rownames(expression_table)) > 1){
      sig_scores[[i]] <- colMeans(expression_table[genes[genes %in% rownames(expression_table)], ])
    }
    if(sum(genes %in% rownames(expression_table)) == 1){
      sig_scores[[i]] <- expression_table[genes[genes %in% rownames(expression_table)], ]
    }
    if(sum(genes %in% rownames(expression_table)) == 0){
      cat(paste0("genes in ", i, " signature not present in expression table\n"))
    }
    rm(i)
  }
  sig_scores <- data.frame(bind_cols(sig_scores))
  rownames(sig_scores) <- colnames(expression_table)
  return(sig_scores)
}
