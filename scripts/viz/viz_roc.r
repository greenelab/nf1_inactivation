# Gregory Way 2016
# NF1 Inactivation Classifier for Glioblastoma
# scripts/viz/viz_roc.r
#
# Visualizes classifier performance by ROC
#
# Usage:
# Run with nf1_roc.py
#
# Output:
# average training/testing ROC and empirical cdf in multiplot

# Load library
library(ggplot2)
library(dplyr)
library(gridExtra)

# Load Command Args
args <- commandArgs(trailingOnly = T)
#args <- c('results/roc_output.tsv', 'figures/GBM_elasticnet_roc_outputalpha_0.15_l1ratio_0.1_roc.pdf')
roc_fh <- args[1]
roc_figure <- args[2]
roc_results <- file.path('results',
                         paste0(tools::file_path_sans_ext(basename(args[2])),
                                '_auroc.tsv'))

# Load Data
roc_data <- readr::read_tsv(roc_fh)

# Save train/test AUROC results
roc_summary <- roc_data %>%
  group_by(type, seed, fold) %>%
  summarise(avg = mean(auc)) %>%
  group_by(type)

mean_auroc <- roc_summary %>% summarise(mean_auroc = mean(avg))
range_auroc_low <- roc_summary %>% summarise(ci_low = quantile(avg, 0.05))
range_auroc_high <- roc_summary %>% summarise(ci_high = quantile(avg, 0.95))

auroc_results <- cbind(mean_auroc, range_auroc_low$ci_low,
                       range_auroc_high$ci_high)
colnames(auroc_results) <- c('type', 'mean', '0.05 CI', '0.95 CI')
write.table(auroc_results, roc_results, row.names = F, sep = '\t')

# Create new variable that stores an ID for each unique iteration
roc_data = roc_data %>% mutate(iteration = paste(seed, fold, sep = '_'))

# Plot
base_theme <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line.x = element_line(size = rel(0.8),
                                               linetype = "solid",
                                               color = "black"),
                    axis.line.y = element_line(size = rel(0.8),
                                               linetype = "solid",
                                               color = "black"),
                    axis.text = element_text(size = rel(0.6), color = "black"),
                    axis.ticks = element_line(size = rel(0.8), color = "black"),
                    axis.ticks.length = unit(0.2, "cm"),
                    axis.title = element_text(size = rel(0.8)),
                    axis.title.y = element_text(vjust = 4.5),
                    axis.title.x = element_text(vjust = -4.5),
                    legend.key.size = unit(0.5, "cm"),
                    legend.text = element_text(size = rel(0.9)),
                    legend.title = element_blank())

roc_grob <- ggplot(roc_data, aes(x = fpr, y = tpr, color = type, fill = type,
                                 group = interaction(seed, fold, type))) +
  labs(x = 'False Positive Rate', y = 'True Positive Rate') + 
  geom_line(stat = "smooth", method = "loess", se = TRUE, alpha = 0.2) + 
  geom_abline(intercept = 0, linetype = "dotted", lwd = rel(0.8)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     limits = c(0, 1.1)) +
  base_theme + theme(plot.margin = unit(c(0.5, 1, 0, 1), 'cm'))

# Summarize the AUROC for each iteration grouped by train/test
roc_auc_data <- roc_data %>% group_by(iteration, type) %>%
  summarize(auc_mean = mean(auc))

dens_cdf_grob <- ggplot(roc_auc_data, aes(x = auc_mean, colour = type,
                                          fill = type)) +
  stat_ecdf(aes(ymin = 0, ymax = ..y..), geom = 'ribbon', alpha = 0.5) +
  labs(x = 'AUROC', y = 'CDF') +
  geom_vline(xintercept = 0.5, linetype = "dotted", lwd = rel(0.8)) +
  base_theme + theme(plot.margin = unit(c(0, 1, 2, 1), 'cm'))

# Extract the legend from the CDF plot
gtable <- ggplot_gtable(ggplot_build(dens_cdf_grob))
legend_grob <- which(sapply(gtable$grobs, function(x) x$name == 'guide-box'))
legend_grob <- gtable$grobs[[legend_grob]]

layout <- matrix(c(rep(1, 90), rep(2, 10)), nrow = 5, ncol = 20)
pdf(roc_figure, height = 7, width = 7)
grid.arrange(arrangeGrob(roc_grob + theme(legend.position = "none"),
                         dens_cdf_grob + theme(legend.position = "none")), 
             layout_matrix = layout,
             legend_grob, nrow = 1)
dev.off()
