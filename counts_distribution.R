################################################################################
#
# counts distribution QC
#
################################################################################

# input is counts matrix: rows is gene names, columns is sample names

### counts distribution --------------------------------------------------------
counts_dataframe_ <- data.frame()
for(i in colnames(counts_distribution_input_)){
  sample_data_ <- data.frame(sample = rep(i, nrow(counts_distribution_input_)),
                             counts = counts_distribution_input_[, i],
                             row.names = NULL)
  counts_dataframe_ <- rbind(counts_dataframe_, sample_data_)
  rm(i, sample_data_)
}
return(counts_dataframe_)
