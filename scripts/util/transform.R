# Gregory Way 2016
# NF1 Inactivation Classifier for Glioblastoma
# scripts/util/transform.R
#
# Transforms the RNAseq distribution to appear more like microarray
#
# Usage:
# Run by transform_rnaseq.py within run_pipeline.sh to transform distribution.
# The script is run with the following arguments:
#
#     R --no-save --args <microarray file> <RNAseq file> <output file>
#
# Output:
# Gene expression matrix for validation set

# Load libraries
library(data.table)

# Load Command Args
args <- commandArgs(trailingOnly = T)

microarray_file <- args[1]
rnaseq_file <- args[2]
output_file <- args[3]

# Load Data
microarray <- data.frame(readr::read_tsv(microarray_file))
rownames(microarray) <- microarray[, 1]
microarray <- microarray[, -1]
rnaseq <- data.frame(readr::read_tsv(rnaseq_file))
rownames(rnaseq) <- rnaseq[, 1]
rnaseq <- rnaseq[, -1]

# Subset microarray to same genes
common_genes <- intersect(rownames(rnaseq), rownames(microarray))
microarray <- microarray[common_genes, ]
rnaseq <- rnaseq[common_genes, ]

# TDM
tcga_tdm <- TDM::tdm_transform(ref_data = 
                                 data.table(cbind(gene = rownames(microarray),
                                                  microarray)),
                               target_data =
                                 data.table(cbind(gene = rownames(rnaseq),
                                                  rnaseq)))
tcga_tdm <- as.matrix(tcga_tdm)
colnames(tcga_tdm)[1] <- "Sample"

# Write the TDM results to file
colnames(tcga_tdm) <- gsub(x = colnames(tcga_tdm), pattern = "[.]",
                           replacement = "-")
write.table(tcga_tdm, output_file, sep = "\t", row.names = F)
