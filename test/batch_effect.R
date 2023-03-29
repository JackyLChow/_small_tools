counts <- readRDS("~/Documents/BFX_proj/_data_public/Chow_PNAS_2020/Chow_PNAS_lognormcounts.rds")
meta <- readRDS("~/Documents/BFX_proj/_data_public/Chow_PNAS_2020/Chow_PNAS_meta_med.rds")

counts_batched <- counts
counts_batched[, 3:8] <- counts_batched[, 3:8] + 15
batch <- factor(c("A", "A", rep("B", 6), rep("A", 8)))

# normal
design_ <- model.matrix(~ meta$treatment)
l_ <- lmFit(counts, design_)
e_ <- eBayes(l_)
deg_1 <- topTable(e_, n = Inf)

# artificial batch w/ adjustment to model
design_ <- model.matrix(~ meta$treatment + batch)
l_ <- lmFit(counts_batched, design_)
e_ <- eBayes(l_)
deg_2 <- topTable(e_, n = Inf)

# artificial batch w/ batch correction
design_ <- model.matrix(~ meta$treatment)
counts_corrected <- removeBatchEffect(counts_batched, batch = batch, design = design_)
l_ <- lmFit(counts_corrected, design_)
e_ <- eBayes(l_)
deg_3 <- topTable(e_, n = Inf)





counts_batched <- counts
counts_batched[, 3:8] <- counts_batched[, 3:8] + 15
counts_batched[, 9:10] <- counts_batched[, 9:10] + 20

batch <- factor(c("A", "A", rep("B", 6), "C", "C", rep("A", 6)))

boxplot(counts_batched)

design_ <- model.matrix(~ meta$treatment)
l_ <- lmFit(counts, design_)
e_ <- eBayes(l_)
deg_1 <- topTable(e_, n = Inf)

design_ <- model.matrix(~ meta$treatment)
l_ <- lmFit(counts_batched, design_)
e_ <- eBayes(l_)
deg_2 <- topTable(e_, n = Inf)

design_ <- model.matrix(~ meta$treatment + batch)
l_ <- lmFit(counts_batched, design_)
e_ <- eBayes(l_)
deg_3 <- topTable(e_, n = Inf)