# Gregory Way 2016
# NF1 Inactivation Classifier for Glioblastoma
# scripts/preprocess_validation.r
#
# Run Quality Control on Validation Microarray and output
# RMA normalized gene expression matrix
#
# Usage:
# Run in run_pipeline.sh and before analyzing validation data
#
# Output:
# Normalized microarray matrix for validation samples

set.seed(123)

#########################################
# Load Libraries
#########################################
library(pd.hta.2.0)
library(hta20stprobeset.db)
library(hta20sttranscriptcluster.db)
library(annotate)
library(affxparser)
library(plyr)

# Every plot in this analysis will be saved to the same pdf file
pdf(file.path('data', 'validation', 'quality_control_assessment.pdf'))

#########################################
# Load and Inspect Data
#########################################
# Read in Cel files
base_dir <- file.path('data', 'validation', 'raw')
celFiles <- c(list.celfiles(file.path(base_dir, 'CelA'), full.names = T),
              list.celfiles(file.path(base_dir, 'CelB'), full.names = T))
affyExpressionFS <- read.celfiles(celFiles)

sample_names <- paste(sapply(sampleNames(affyExpressionFS),
                             function(x) {substr(x, 4, 6)}))

# Visually inspect for spatial artifacts
for (samp in 1:ncol(affyExpressionFS)) {
  oligo::image(affyExpressionFS[, samp], transfo = rank)
}

# MA Plots
oligo::MAplot(affyExpressionFS)

#########################################
# Process Data
#########################################
# RMA correction
exp.rma <- oligo::rma(affyExpressionFS, target = "core")
exp.rma_probe <- oligo::rma(affyExpressionFS, target = "probeset")

# MA Plots after RMA
oligo::MAplot(exp.rma)

# Box Plots
oligo::boxplot(exp.rma, names = sample_names)

# Inspect distribution of data
oligo::hist(exp.rma, transfo = identity)

# Get NetAffx Biological Annotation and store it in the feature slot
featureData(exp.rma) <- getNetAffx(exp.rma, "transcript")
featureData(exp.rma_probe) <- getNetAffx(exp.rma_probe, "probeset")

# Get Probes
exp.trans <- exprs(exp.rma)
exp.probe <- exprs(exp.rma_probe)

# Use HTA2.0 Transcript probes for downstream analyses

#########################################
# Is Batch Correction Necessary?
# Use Guided PCA to determine: http://www.ncbi.nlm.nih.gov/pubmed/23958724
#########################################
batch <- c(rep(1, 6), rep(2, 6))  # For plate number
batch_info <- c()
batch_date <- c()

# For sample processing date
for (cel in celFiles) {
  binfo <- readCelHeader(cel)$datheader
  batch_info <- rbind(batch_info, binfo)
  bdate <- gsub(".*([0-9]{2,2}/[0-9]{2,2}/[0-9]{2,2}).*","\\1",binfo)
  batch_date <- c(batch_date, bdate)
}

batch_date_ready <- batch_date
for (i in 1:length(unique(batch_date))) {
  batch_date_ready[batch_date == unique(batch_date)[i]] <- as.integer(i)
}

# Guided PCA to detect batch
batch_out <- gPCA::gPCA.batchdetect(x = t(exp.trans), batch, scaleY = T)
batch_out_date <- gPCA::gPCA.batchdetect(x = t(exp.trans),
                                         as.integer(batch_date_ready),
                                         scaleY = T)

print('Does Batch Have a significant effect: Plate')
batch_out$delta  # [1] 0.2197712
batch_out$p.val  # [1] 0.744

print('Does Batch Have a significant effect: Date Processed')
batch_out_date$delta  # [1] 0.3493635
batch_out_date$p.val  # [1] 0.79

# Plot Results
gPCA::PCplot(batch_out, ug = "unguided", type = "1v2")
gPCA::PCplot(batch_out_date, ug = "unguided", type = "1v2")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Based on this analysis, there is no need to Batch Correct
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#########################################
# Extract info to map to genes
#########################################
# Get Probe Mappings
x <- hta20stprobesetSYMBOL

# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- unlist(as.list(x[mappedkeys(x)]))

# Get Transcript Mappings
x <- hta20sttranscriptclusterSYMBOL

# Get the probe identifiers that are mapped to a gene symbol
mapped_transcripts <- unlist(as.list(x[mappedkeys(x)]))
TranscriptLevel <- exp.trans[names(mapped_transcripts), ]

# Process gene level data data
Gene <- cbind(mapped_transcripts, TranscriptLevel)
Gene <- Gene[order(Gene[, 1]), ]
Gene <- cbind(Gene[, 1], rownames(Gene), Gene[, 2:ncol(Gene)])

# To map transcripts to gene, take the maxmean from WGCNA collapseRows function
output_data <- WGCNA::collapseRows(TranscriptLevel,
                                   rowGroup = paste(mapped_transcripts),
                                   rowID = names(mapped_transcripts))

# Write validation set data
val_file = file.path('data', 'validation', 'normalized', 'validation_set.tsv')
write.table(output_data$datETcollapsed, val_file, sep = "\t", row.names = T,
            col.names = NA)
dev.off()
