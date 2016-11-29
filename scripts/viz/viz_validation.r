# Gregory Way 2016
# NF1 Classifier GBM/LGG
# scripts/viz/viz_validation.R
#
# vizualize the validation results
#
# Usage:
# Run by run_pipeline.sh
#
# Output:
# Four validation set figures and one table
# 1) Boxplot of ensemble predictions
# 2) Scatterplot of predictions by quantified protein
# 3) Boxplot of predictions for NF1 WT/Inactive protein
# 4) Power analysis plot of n x power
# 5) Results of t-test/effect size analysis

# Load libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(calibrate)
library(pwr)
library(lsr)

# Load Command Args
args <- commandArgs(trailingOnly = T)
print(args)
validation_fh <- args[1]
protein_fh <- args[2]
figure_main <- paste0(args[3], '_main_figure.pdf')
figure_power <- paste0(args[3], '_power_analysis.pdf')
table_stat <- paste0(args[3], '_ttest_power_table.csv')

# Load Constants and Data
qc_pass <- c('CB2', 'H5M', '3HQ', 'PBH', 'LNA', 'YXL', 'VVN', 'R7K',
             'RIW', 'TRM', 'UNY', 'W31')
validation <- readr::read_tsv(validation_fh)
protein <- read.csv(protein_fh, header = T, stringsAsFactors = F)

# Rename samples
validation$sample_id <- apply(validation, 1, function(x) {substr(x[1], 4, 6)})

# Subset validation set to only quality samples
validation <- validation[validation$sample_id %in% qc_pass, ]

# Make 0 estimates -1
validation$prediction[validation$prediction == 0] <- -1

# Finalize ensemble predictions (used later for final call)
sample_groups <- validation %>% group_by(sample_id)
neg_prob <- sample_groups %>% summarise(neg_prob = sum(neg_prob))
pos_prob <- sample_groups %>% summarise(pos_prob = sum(pos_prob))
pred_count <- sample_groups %>% summarise(pred_sum = sum(prediction))

# Create new column of predictions weighted by test set AUROC
validation <- validation %>%  mutate(weighted_auc = prediction * test_auc)
weight_auc <- validation %>% group_by(sample_id) %>% 
  summarise(weight_auc = mean(weighted_auc))

# Obtain protein by prediction correlations
plot_ready <- cbind(pred_count, pos_prob, neg_prob, weight_auc,
                    protein[match(pred_count$sample_id, protein$sample_id), ])
plot_ready <- plot_ready[complete.cases(plot_ready), ]
plot_ready <- plot_ready[, c('sample_id', 'pred_sum', 'pos_prob', 'neg_prob',
                             'weight_auc', 'u87pi.norm', 'plate')]

# Relabel to add plate info
plot_ready[13:14, 6:7] <- protein[12:13, c(2, 4)]
plot_ready[13:14, 1:5] <- plot_ready[c(2, 3), 1:5]
colnames(plot_ready)[ncol(plot_ready)] <- 'Plate'

# Get the mean values for samples measured on both plates
h5m_mean <- mean(plot_ready[plot_ready$sample_id == 'H5M', "u87pi.norm"])
cb2_mean <- mean(plot_ready[plot_ready$sample_id == 'CB2', "u87pi.norm"])

# Subset data to find correlations
protein_values <- plot_ready[1:(nrow(plot_ready) - 2), ]
protein_values[protein_values$sample_id == 'H5M', "u87pi.norm"] <- h5m_mean
protein_values[protein_values$sample_id == 'CB2', "u87pi.norm"] <- cb2_mean
protein_values <- protein_values %>%
  mutate(protein_rank = dense_rank(desc(u87pi.norm)))

# Refactor and match protein concentration for plotting
validation$sample_id <- factor(validation$sample_id, levels = qc_pass)
validation <- validation %>% rowwise() %>% 
  mutate(protein = protein_values[protein_values$sample_id == sample_id, 
                                  'protein_rank'])

# Plot Paddle Chart
paddle_grob <- ggplot(validation, aes(x = sample_id, y = weighted_auc,
                                      fill = protein)) +
  geom_violin(adjust = 0.8) +
  xlab('') + ylab('Weighted Prediction Scores') +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.5) +
  #scale_fill_brewer(palette = 'Greys') +
  scale_fill_distiller(palette = "Blues") +
  theme(title = element_text(size = rel(2.2)),
        axis.title = element_text(size = rel(0.4)),
        axis.text.x = element_text(size = rel(1), angle = 45),
        axis.text.y = element_text(size = rel(1)),
        legend.position = "none",
        panel.grid.major = element_line(color = 'white', size = 0.3),
        panel.grid.minor = element_line(color = 'white', size = 0.3),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = 'black', size = 0.5),
        axis.ticks = element_blank())

# Plot scatter plot
scatter_grob <- ggplot(plot_ready, aes(x = weight_auc, y = u87pi.norm,
                                       color = Plate, label = sample_id)) +
  xlab('Weighted Predictions (Mean)') +
  ylab('NF1 Protein Concentration') +
  geom_point(size = rel(0.8)) +
  geom_text(color = 'black', hjust = 1.2, size = rel(2.2)) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.5), limits = c(-1, 1)) +
  theme(title = element_text(size = rel(2.2)),
        axis.title = element_text(size = rel(0.4)),
        axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.line.x = element_line(color = 'black', size = 0.5),
        axis.line.y = element_line(color = 'black', size = 0.5),
        legend.title = element_text(size = rel(0.5)),
        legend.key = element_rect(fill = 'white'),
        legend.text = element_text(size = rel(0.6)),
        legend.margin = unit(0.05, "cm"), 
        panel.grid.major = element_line(color = 'white', size = 0.3),
        panel.grid.minor = element_line(color = 'white', size = 0.3),
        panel.background = element_rect(fill = 'white'))

# Plot barchart for positive and negative predictions
predict_nf1 <- vector(length = nrow(plot_ready))
predict_nf1[plot_ready$weight_auc > 0] <- 'NF1 Inactive'
predict_nf1[plot_ready$weight_auc <= 0] <- 'NF1 Wildtype'

# Get prediction plot ready
predict_plot <- cbind(plot_ready, predict_nf1)

# Refactor for plotting
predict_plot$predict_nf1 <- factor(predict_plot$predict_nf1,
                                   levels = c('NF1 Wildtype', 'NF1 Inactive'))

# Reassign CB2 and H5M
predict_plot <- predict_plot[1:(nrow(predict_plot) - 2), ]
predict_plot[predict_plot$sample_id == 'H5M', 'u87pi.norm'] <- h5m_mean
predict_plot[predict_plot$sample_id == 'CB2', 'u87pi.norm'] <- cb2_mean

# Set seed for geom_jitter()
set.seed(1234)
box_plot_grob <- ggplot(predict_plot, aes(x = predict_nf1, y = u87pi.norm)) +
  geom_boxplot(aes(fill = predict_nf1), lwd = 0.2, outlier.colour = 'white') +
  xlab('') + ylab('NF1 Protein Concentration') +
  geom_jitter(width = 0.4, size = rel(0.8)) +
  scale_fill_manual(values = c('lightblue', 'red')) +
  theme(title = element_text(size = rel(2.2)),
        axis.title = element_text(size = rel(0.4)),
        axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        legend.position = "none",
        panel.grid.major = element_line(color = 'white', size = 0.3),
        panel.grid.minor = element_line(color = 'white', size = 0.3),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color = 'black', size = 0.5),
        axis.line.y = element_line(color = 'black', size = 0.5),
        axis.ticks.x = element_blank())

# Arrage all grobs and save figure
blankGrob <- ggplot() + geom_blank() + theme(panel.background = element_blank())
layout <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
pdf(figure_main, height = 5.5, width = 6.5)
grid.arrange(blankGrob, paddle_grob, scatter_grob, box_plot_grob,
             layout_matrix = layout)
dev.off()

# t-test for difference in protein between groups
nf1_inactive <- predict_plot[predict_plot$predict_nf1 == 'NF1 Inactive',
                             'u87pi.norm']
nf1_wildtype <- predict_plot[predict_plot$predict_nf1 == 'NF1 Wildtype',
                             'u87pi.norm']
results <- t.test(x = nf1_inactive, y = nf1_wildtype, alternative = 'less')

# Store results in table
t_stat <- results$statistic
p_val <- results$p.value
effect_size <- cohensD(nf1_inactive, nf1_wildtype)
required_n <- round(pwr.t.test(d = effect_size, power = 0.8 ,
                               type = 'one.sample')$n)
full_results <- cbind(t_stat, p_val, effect_size, required_n)
write.table(full_results, table_stat, sep = ',', row.names = F)

# Power analysis figure
n <- seq(2, 100)
power <- sapply(n, function(x) pwr.t.test(d = effect_size, n = x,
                                          type = 'one.sample')$power)
power_plot <- data.frame(cbind(n, power))
pdf(figure_power, height = 3, width = 3)
ggplot(data = power_plot, aes(x = n, y = power)) + geom_line(lwd = 0.5) +
  theme(axis.title = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(0.8)),
        axis.text.y = element_text(size = rel(0.8)),
        legend.position = "none",
        panel.grid.major = element_line(color = 'white', size = 0.3),
        panel.grid.minor = element_line(color = 'white', size = 0.3),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color = 'black', size = 0.5),
        axis.line.y = element_line(color = 'black', size = 0.5)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", lwd = 0.5, col = 'red') +
  geom_vline(xintercept = 12, linetype = 'dashed', lwd = 0.5, col = 'blue')
dev.off()