# (C) Gregory Way 2016
# NF1 Classifier GBM/LGG
# scripts/viz/viz_cv.r
#
# Visualizes cross validation results
#
# Usage:
# Run with nf1_classifier.py
#
# Output:
# Training and testing error plots for hyperparameter selection

# Load library
library(ggplot2)
library(dplyr)

# Load Command Args
args <- commandArgs(trailingOnly = T)
print(args)
data_fh <- args[1]
output_fh <- args[2]
ml <- args[3]

# Load Data
data_cv <- readr::read_tsv(data_fh)

# Multiplot name
multiplot_fh <- paste0(unlist(strsplit(output_fh, '[.]'))[1], '_multiplot.png')

# Construct classifier specific plots
if (ml == 'SVM') {
  x_label <- 'Hyperparameter: C'
  param_sweep <- ggplot(data_cv, aes(x = param, y = error, colour = class, fill = class))
  h <- 400
  w <- 600
  title <- ''
} else {
  x_label <- 'Hyperparamter: l1 Ratio'
  w <- 700
  h <- 500
  title <- 'Hyperparameter: Alpha'
  param_sweep <- ggplot(data_cv, aes(x = l1_ratio, y = error, colour = class, fill = class)) + 
    facet_wrap(~ alpha)
}

##############################
# Plot Parameter Sweep
##############################
param_sweep <- param_sweep + stat_smooth(method = 'loess', se = TRUE) + 
  labs(x = x_label, y = 'AUROC') +
  ggtitle(title) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        title = element_text(size = rel(2)),
        axis.line.x = element_line(size = rel(1), linetype = "solid", color = "black"), 
        axis.line.y = element_line(size = rel(1), linetype = "solid", color = "black"), 
        axis.text = element_text(size = rel(1), color = "black"),
        axis.ticks = element_line(size = rel(1), color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = rel(1.2)),
        axis.title.y = element_text(vjust = 4.5),
        axis.title.x = element_text(vjust = -4.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)), 
        legend.title = element_text(size = rel(1.2)))

# Save
png(output_fh, height = h, width = w)
param_sweep
dev.off()

##############################
# Data subset (Plot Parameter Selection)
##############################
# What is best parameters for test set?
if (ml == 'SVM') {
  param <- data_cv[data_cv$class == 'test', ] %>% group_by(alpha) %>% summarise(avg = mean(error))
} else if (ml == 'elasticnet') {
  param <- data_cv[data_cv$class == 'test', ] %>% group_by(alpha, l1_ratio) %>% summarise(avg = mean(error))
}

# What are the best parameters with the lowest test set accuracy?
best_parameters <- param[param$avg == max(param$avg), ]

# Plot Best Alpha Range
best_alpha <- best_parameters$alpha
data_sub <- data_cv[data_cv$alpha == best_alpha, ]

best_alpha_plot <- ggplot(data_sub, aes(class, error)) +
  stat_boxplot(aes(fill = class)) +
  xlab('') + ylab('AUROC') + facet_grid(. ~ l1_ratio) +
  theme(title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        legend.key = element_rect(fill = 'white'),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        panel.grid.major = element_line(color = 'grey', size = 0.3),
        panel.grid.minor = element_line(color = 'grey', size = 0.3),
        panel.background = element_rect(fill = 'white'),
        axis.line = element_line(color = 'black', size = 0.5),
        axis.ticks = element_blank(),
        title = element_text(size = 25)) +
       # plot.margin = unit(c(1.2, 2, 0, 1.5), 'cm')) + 
  ggtitle(paste0('l1 Ratio', '\nalpha = ', best_alpha))

# Save
output_fh <- paste0(unlist(strsplit(output_fh, '[.]'))[1], '_alpha_', best_alpha, '.png')
png(output_fh, height = h, width = w)
best_alpha_plot
dev.off()

##############################
# Plot best overall parameters
##############################
# Subset to best l1_ratio if elasticnet
best_l1 <- best_parameters$l1_ratio

l1_sub <- data_sub[data_sub$l1_ratio == best_l1, ]

best_l1_plot <- ggplot(l1_sub, aes(class, error)) +
  stat_boxplot(aes(fill = class)) +
  xlab('') + ylab('AUROC') +
  theme(title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        legend.key = element_rect(fill = 'white'),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        panel.grid.major = element_line(color = 'grey', size = 0.3),
        panel.grid.minor = element_line(color = 'grey', size = 0.3),
        panel.background = element_rect(fill = 'white'),
        axis.line = element_line(color = 'black', size = 0.5),
        axis.ticks = element_blank(),
        title = element_text(size = 25)) +
       # plot.margin = unit(c(1.2, 2, 0, 1.5), 'cm')) +
  ggtitle(paste0('l1 ratio = ', best_l1, '\nalpha = ', best_alpha))

# Save
output_fh <- paste0(unlist(strsplit(output_fh, '[.]'))[1], best_alpha,
                    '_l1ratio_', best_l1, '.png')
png(output_fh, height = 270, width = 300)
best_l1_plot
dev.off()

##############################
# Save Multiplot
##############################
all_plots <- list(param_sweep + theme(legend.position = "none"), best_alpha_plot, best_l1_plot)
layout <- matrix(c(rep(1, 20), rep(c(2, 2, 2, 3, 3), 4)), nrow = 5, ncol = 8)

png(multiplot_fh, height = 500, width = 950)
Rmisc::multiplot(plotlist = all_plots, layout = layout)
dev.off()