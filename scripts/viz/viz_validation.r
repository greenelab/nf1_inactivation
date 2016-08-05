# (C) Gregory Way 2016
# NF1 Classifier GBM/LGG
# scripts/viz/viz_validation.r
#
# vizualize the validation results
#
# Usage:
# Run by ANALYSIS.sh
#
# Output:
# Three validation set figures
# 1) Boxplot of ensemble predictions
# 2) Scatterplot of predictions by quantified protein
# 3) Boxplot of predictions for NF1 WT/Inactive protein

# Load libraries
library(ggplot2)
library(dplyr)
library(calibrate)

# Load Command Args
args <- commandArgs(trailingOnly = T)
print(args)
validation_fh <- args[1]
protein_fh <- args[2]
figure_box <- paste0(args[3], 'prediction_box.png')
figure_scatter <- paste0(args[3], 'protein_scatter.png')
figure_pred <- paste0(args[3], 'protein_box.png')

# Load Constants and Data
qc_pass <- c('CB2', 'H5M', '3HQ', 'PBH', 'LNA', 'YXL', 'VVN', 'R7K', 'RIW', 'TRM', 'UNY', 'W31')
validation <- readr::read_tsv(validation_fh)
protein <- read.csv(protein_fh, header = T, stringsAsFactors = F)

# Rename samples
validation$sample_id <- apply(validation, 1, function(x) {substr(x[1], 4, 6)})

# Subset validation set to only quality samples
validation <- validation[validation$sample_id %in% qc_pass, ]

# Make 0 estimates -1
validation$prediction[validation$prediction == 0] <- -1

# Finalize ensemble predictions (used later for final call)
neg_prob <- validation %>% group_by(sample_id) %>% summarise(neg_prob = sum(neg_prob))
pos_prob <- validation %>% group_by(sample_id) %>% summarise(pos_prob = sum(pos_prob))
pred_count <- validation %>% group_by(sample_id) %>% summarise(pred_sum = sum(prediction))

# Create new column of predictions weighted by test set AUROC
validation <- validation %>%  mutate(weighted_auc = prediction * test_auc)
weight_auc <- validation %>% group_by(sample_id) %>% summarise(weight_auc = mean(weighted_auc))

# Refactor for plotting
validation$sample_id <- factor(validation$sample_id, levels = qc_pass)

###################################
# Plot Bar Chart
###################################
png(figure_box, height = 250, width = 360)
ggplot(validation, aes(x = sample_id, y = weighted_auc, fill = sample_id)) +
  geom_violin(adjust = 0.8) + xlab('') + ylab('Weighted AUROC') +
  geom_hline(yintercept = 0, linetype="dashed", lwd = 1) +
  theme(title = element_text(size = rel(2.2)),
        axis.title = element_text(size = rel(0.6)),
        axis.text.x = element_text(size = rel(1.5), angle = 45),
        axis.text.y = element_text(size = rel(1.5)),
        legend.position = "none",
        panel.grid.major = element_line(color = 'white', size = 0.3),
        panel.grid.minor = element_line(color = 'white', size = 0.3),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = 'black', size = 0.5),
        axis.ticks = element_blank())
dev.off()

###################################
# Plot scatter plot
###################################
plot_ready <- cbind(pred_count, pos_prob, neg_prob, weight_auc,
                    protein[match(pred_count$sample_id, protein$X), ])
plot_ready <- plot_ready[complete.cases(plot_ready), ]
plot_ready <- plot_ready[, c('sample_id', 'pred_sum', 'pos_prob', 'neg_prob', 'weight_auc',
                             'u87pi.norm', 'plate')]

# Relabel to add plate info
plot_ready[13:14, 6:7] <- protein[12:13, c(2, 4)]
plot_ready[13:14, 1:5] <- plot_ready[c(2, 3), 1:5]
colnames(plot_ready)[ncol(plot_ready)] <- 'Plate'

# Get the mean values for samples measured on both plates
h5m_mean <- mean(plot_ready[plot_ready$sample_id == 'H5M', "u87pi.norm"])
cb2_mean <- mean(plot_ready[plot_ready$sample_id == 'CB2', "u87pi.norm"])

# Subset data to find correlations
cor_test_ready <- plot_ready[1:(nrow(plot_ready) - 2), ]
cor_test_ready[cor_test_ready$sample_id == 'H5M', "u87pi.norm"] <- h5m_mean
cor_test_ready[cor_test_ready$sample_id == 'CB2', "u87pi.norm"] <- cb2_mean

# Spearman correlation
results <- cor.test(x = cor_test_ready$weight_auc, y = cor_test_ready$u87pi.norm,
                    method = 'spearman')

# Plot
add_text <- paste0('rho = ', round(results$estimate, 2) , '\np = ', round(results$p.value, 2))
png(figure_scatter, height = 250, width = 360)
ggplot(plot_ready, aes(x = weight_auc, y = u87pi.norm, color = Plate, label = sample_id)) +
  xlab('Weighted Predictions (Mean)') + ylab('NF1 Protein Level') +
  geom_point() + geom_text(color = 'black', vjust = -0.3, size = rel(3)) +
  theme(title = element_text(size = rel(2)),
        axis.title = element_text(size = rel(0.6)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(0.6)),
        legend.key = element_rect(fill = 'white'),
        legend.text = element_text(size = rel(1)),
        panel.grid.major = element_line(color = 'white', size = 0.3),
        panel.grid.minor = element_line(color = 'white', size = 0.3),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color = 'black', size = 0.5),
        axis.line.y = element_line(color = 'black', size = 0.5),
        axis.ticks = element_blank(),
        title = element_text(size = 25))
dev.off()

###################################
# Plot barchart for positive and negative predictions
###################################
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

png(figure_pred, height = 250, width = 360)
ggplot(predict_plot, aes(x = predict_nf1, y = u87pi.norm)) +
  geom_boxplot(aes(fill = predict_nf1), outlier.colour = 'white') +
  xlab('') + ylab('NF1 Protein Level') +
  geom_jitter(width = 0.2) +
  scale_fill_manual(values = c('lightblue', 'red')) +
  theme(title = element_text(size = 20),
        axis.title = element_text(size = rel(0.6)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        legend.position = "none",
        panel.grid.major = element_line(color = 'white', size = 0.3),
        panel.grid.minor = element_line(color = 'white', size = 0.3),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color = 'black', size = 0.5),
        axis.line.y = element_line(color = 'black', size = 0.5),
        axis.ticks = element_blank(),
        title = element_text(size = 25))
dev.off()
