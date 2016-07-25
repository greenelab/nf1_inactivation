# (C) Gregory Way 2016
# NF1 Classifier GBM/LGG
# scripts/viz/viz_geneweights.r
#
# vizualize the distribution of gene predictive weights
#
# Usage:
# Run by ANALYSIS.sh 
#
# Output:
# 1) Figure - distribution of gene weights across ensemble
# 2) Tables - high weight genes (positive and negative contributions)

set.seed(123)

# Load libraries
library(ggplot2)
library(dplyr)
library(calibrate)

# Load Command Args
args <- commandArgs(trailingOnly = T)
print(args)
weight_fh <- args[1]
out_fh <- args[2]

# Load Constants and Data
weight <- read.table(weight_fh, header = F, stringsAsFactors = F)
colnames(weight) <- c('gene', 'weight')
weight <- weight[order(weight$weight, decreasing = F), ]
weight$rank <- 1:nrow(weight)

# Scale the weights to unit norm
weight$weight <- weight$weight / nrow(weight)

# Define high weight genes by number of standard deviations
high_weight <- sd(weight$weight) * 2

# Plot
png(out_fh, height = 370, width = 550)
ggplot(weight, aes(x = rank, y = weight)) + 
  geom_point() + xlab('Rank') + ylab('Gene Weight') +
  geom_hline(yintercept = 0, linetype = "dotted", lwd = 1, color = 'red') +
  geom_hline(yintercept = high_weight, linetype = "dashed", lwd = 1.2, color = 'red') +
  geom_hline(yintercept = -high_weight, linetype = "dashed", lwd = 1.2, color = 'red') +
  theme(title = element_text(size = rel(2.2)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        legend.position = "none",
        panel.grid.major = element_line(color = 'white', size = 0.3),
        panel.grid.minor = element_line(color = 'white', size = 0.3),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color = 'black', size = 0.5),
        axis.line.y = element_line(color = 'black', size = 0.5))
dev.off()

# Output text files of high weight positive and high weight negative genes
high_weight_pos <- weight[weight$weight > high_weight, ]
high_weight_pos <- high_weight_pos[order(high_weight_pos$weight, decreasing = T), ]
high_weight_pos$rank <- seq(1, nrow(high_weight_pos))

high_weight_neg <- weight[weight$weight < -high_weight, ]
high_weight_neg$rank <- seq(1, nrow(high_weight_neg))

write.table(high_weight_pos, 'results/high_weight_genes_positive.tsv', sep = '\t', row.names = F)
write.table(high_weight_neg, 'results/high_weight_genes_negative.tsv', sep = '\t', row.names = F)
