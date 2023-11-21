# fold_changer

# input is two tables:
## 1) matrix, table_of_values to be converted to fold change vs. reference
### rows is observations, columns is sample_id
## 2) data.frame, table_of_samples
### columns are: sample_id, subject, visit
## 3) character, reference_visit

fold_changer <- function(table_of_values = l2tpm[c("CD8A", "CD8B"), ],
                         table_of_samples = metadata,
                         sample_column = "barcode",
                         subject = "subject",
                         visit = "visit",
                         reference_visit = "Screening"){
  reference_samples <- table_of_samples[table_of_samples[, visit] == reference_visit, ]
  subjects <- unique(reference_samples[, subject])[1:2]
  results <- list()
  for(i in subjects){
    subject_table <- table_of_samples[table_of_samples[, subject] == i, ]
    subject_visits <- subject_table[, visit]
    subject_samples <- subject_table[, sample_column]
    subject_reference_sample <- subject_table[subject_table[, visit] == reference_visit, sample_column]
    
    
    
    results[[i]] <- data.frame(subject = i, visit = subject_visits,
                               reference_sample = subject_reference_sample,
                               sample = subject_samples)
  }
  print(do.call(rbind, results))
}

fold_changer()
