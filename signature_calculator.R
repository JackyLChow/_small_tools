# signature scores
## calculate scores
sig_scores_ <- list()
for(i in names(signatures)){
  sig_s_ <- log2(counts_quantile_normalized[signatures[[i]], metadata$sequencingID] + 1)
  sig_scores_[[i]] <- colMeans(sig_s_)
  rm(sig_s_, i)
}