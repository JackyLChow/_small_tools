# fold_changer

# input is two tables:
## 1) matrix, table_of_values to be converted to fold change vs. reference
### rows is observations, columns is sample_id
## 2) data.frame, table_of_samples
### columns are: sample_id, subject, visit
## 3) character, reference_visit

# try to keep in base R

fold_changer <- function(table_of_values = l2tpm[c("CD8A", "CD8B"), ],
                         table_of_samples = metadata,
                         sample_column = "sample",
                         subject = "subject",
                         visit = "visit",
                         reference_visit = "Screening"){
  # viable subjects must have multiple samples
  viable_subjects <- unique(table_of_samples[duplicated(table_of_samples[, subject]), subject])
  # identify referable samples
  reference_samples <- table_of_samples[table_of_samples[, visit] == reference_visit &
                                          table_of_samples[, subject] %in% viable_subjects, ]
  # filter viable subjects with reference samples
  subjects <- unique(reference_samples[, subject])
  results <- list()
  for(i in subjects){
    print(i)
    # parse out subject level data
    subject_table <- table_of_samples[table_of_samples[, subject] == i, ]
    subject_visits <- subject_table[, visit]
    subject_samples <- subject_table[, sample_column]
    subject_reference_sample <- subject_table[subject_table[, visit] == reference_visit, sample_column]
    
    # extract subject values
    reference_values <- table_of_values[, subject_reference_sample]
    print(dim(reference_values))
    
    # calculate log2 fold changes
    l2fc <- list()
    for(j in subject_samples){
      print(j)
      l2fc[[j]] <- log2(table_of_values[, j] + 1) - log2(reference_values + 1)
    }
    l2fc <- do.call(rbind, l2fc)
    
    # compile final table
    results[[i]] <- data.frame(subject = i, visit = subject_visits,
                               reference_sample = subject_reference_sample,
                               sample = subject_samples,
                               l2fc)
  }
  results <- do.call(rbind, results)
  return(results)
}

foo <- fold_changer()
