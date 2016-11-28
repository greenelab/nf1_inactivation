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
args <- c('results/roc_tdm_output.tsv', 'figures/test_roc_agg_tdm.pdf')
roc_fh <- args[1]
roc_figure <- args[2]
roc_results <- file.path("results",
                         paste0(tools::file_path_sans_ext(basename(args[2])),
                                "_auroc.tsv"))

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
colnames(auroc_results) <- c("type", "mean", "0.05 CI", "0.95 CI")
write.table(auroc_results, roc_results, row.names = F, sep = "\t")

# Create new variable that stores an ID for each unique iteration
roc_data <- roc_data %>% mutate(iteration = paste(seed, fold, sep = "_"))

# Aggregate mean false positive rates for each iteration
roc_data_aggregate <- roc_data %>%
  group_by(seed, fold, type, tpr) %>%
  summarise(mean_fpr = mean(fpr)) %>%
  group_by(seed, fold, type) %>%
  mutate(full_group = paste(seed, fold, type, sep = "_"))

# When ROC curve iterations are saved, they often lack a zero element
# for both true positives and true negatives. Add the point (0, 0) when
# this is the case.
update_rows <- list()
index <- 1
for (group in unique(roc_data_aggregate$full_group)) {
  subset_group <- roc_data_aggregate[roc_data_aggregate$full_group == group, ]
  if (nrow(subset_group[subset_group$mean_fpr == 0 |
                        subset_group$tpr == 0, ]) == 0) {
    add_row <- subset_group[1, ]
    add_row$tpr <- add_row$mean_fpr <- 0
    update_rows[[index]] <- add_row
    index <- index + 1
  }
}

# Add (0, 0) points to the aggregate dataframe
roc_data_aggregate <- dplyr::bind_rows(roc_data_aggregate,
                                       do.call(rbind, update_rows))

# Take the mean of the aggregate by iteration
roc_data_mean <- roc_data_aggregate %>% group_by(type, tpr) %>%
  summarize(full_mean_fpr = mean(mean_fpr))

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

roc_grob <- ggplot(roc_data_mean, aes(x = full_mean_fpr, y = tpr,
                                           color = type, fill = type)) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  geom_line(size = rel(0.6)) + geom_point(size = rel(0.4)) +
  geom_abline(intercept = 0, linetype = "dashed", lwd = rel(0.8)) +
  geom_line(data = roc_data_aggregate, inherit.aes = FALSE,
            aes(x = mean_fpr, y = tpr, color = type, fill = type,
                group = interaction(seed, fold, type)), 
            alpha = 0.1, size = 0.1) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     limits = c(0, 1)) +
  base_theme + theme(plot.margin = unit(c(0.1, 0, 0, 0), "cm"))

# Summarize the AUROC for each iteration grouped by train/test
roc_auc_data <- roc_data %>% group_by(iteration, type) %>%
  summarize(auc_mean = mean(auc))

dens_cdf_grob <- ggplot(roc_auc_data, aes(x = auc_mean, colour = type,
                                          fill = type)) +
  stat_ecdf(aes(ymin = 0, ymax = ..y..), geom = "ribbon", alpha = 0.5) +
  labs(x = "AUROC", y = "CDF") +
  geom_vline(xintercept = 0.5, linetype = "dashed", lwd = rel(0.8)) +
  base_theme + theme(plot.margin = unit(c(0.8, 0, 0, 0), "cm"))

# Extract the legend from the CDF plot
gtable <- ggplot_gtable(ggplot_build(dens_cdf_grob))
legend_grob <- which(sapply(gtable$grobs, function(x) x$name == "guide-box"))
legend_grob <- gtable$grobs[[legend_grob]]

layout <- matrix(nrow = 50, ncol = 50)
layout[1:32, 1:42] <- 1
layout[33:50, 1:42] <- 2
layout[1:50, 43:50] <- 3

pdf(roc_figure, height = 5, width = 4.5)
grid.arrange(roc_grob + theme(legend.position = "none"),
             dens_cdf_grob + theme(legend.position = "none"),
             layout_matrix = layout,
             legend_grob, nrow = 1)
dev.off()
