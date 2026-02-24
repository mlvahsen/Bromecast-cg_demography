# This code creates R2 and confusion matrix plots to assess model accuracy

library(tidyverse); library(geomtextpath); library(cmdstanr)
theme_set(theme_classic(base_size = 16))

# Read in model output
fit <- as_cmdstan_fit(files = c("outputs/demo_model_fecun-202506091642-1-6d200d.csv",
                                "outputs/demo_model_fecun-202506091642-2-6d200d.csv",
                                "outputs/demo_model_fecun-202506091642-3-6d200d.csv"))

fit_s <- as_cmdstan_fit(files = c("outputs/demo_model_surv_noncenter-202506091738-1-51da52.csv",
                                  "outputs/demo_model_surv_noncenter-202506091738-2-51da52.csv",
                                  "outputs/demo_model_surv_noncenter-202506091738-3-51da52.csv"))

# Get model draws
draws_f <- fit$draws(format = "df")
draws_s <- fit_s$draws(format = "df")

# Read in model data
cg_model <- read_csv("data/cg_model_data.csv")

# Make data set of just plants that survived and reproduced
cg_model$survived <- ifelse(cg_model$seed_count > 0, 1, 0)

# Make subset of data that survived
cg_model %>% 
  filter(survived == 1) -> cg_model_fecun

y_rep_fecun <- draws_f %>% select(contains("seed_count_pred")) %>% as.matrix()
ppc_means_fecun <- rowMeans(y_rep_fecun)
ppc_sds_fecun <- apply(y_rep_fecun, 1, sd)
obs_mean_fecun <- mean(cg_model_fecun$seed_count)
obs_sd_fecun <- sd(cg_model_fecun$seed_count)

# PPC for mean
tibble(mean = ppc_means_fecun) %>% 
  ggplot(aes(x = mean)) +
  geom_density(fill = "dodgerblue", alpha = 0.5, color = "darkblue", linewidth = 2) +
  xlab("posterior predictive mean") +
  geom_vline(aes(xintercept = obs_mean_fecun),
                               color = "dodgerblue", linewidth = 2,
                               linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(ppc_means_fecun, 0.025)), color = "darkblue",
             linetype = "dotted") +
  geom_vline(aes(xintercept = quantile(ppc_means_fecun, 0.975)), color = "darkblue",
             linetype = "dotted")-> ppc_mean_plot

# PPC for sd
tibble(sd = ppc_sds_fecun) %>% 
  ggplot(aes(x = sd)) +
  geom_density(fill = "gold", alpha = 0.5, color = "goldenrod", linewidth = 2) +
  xlab("posterior predictive std. dev.") +
  geom_vline(aes(xintercept = obs_sd_fecun),
                               color = "gold", linewidth = 2,
                               linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(ppc_sds_fecun, 0.025)), color = "goldenrod",
             linetype = "dotted") +
  geom_vline(aes(xintercept = quantile(ppc_sds_fecun, 0.975)), color = "goldenrod",
             linetype = "dotted") -> ppc_sd_plot
  
# Survival PPC
y_rep_surv <- draws_s %>% select(contains("survival_pred")) %>% as.matrix()
ppc_success <- rowSums(y_rep_surv)
obs_success <- sum(cg_model$survived)

tibble(success = ppc_success) %>% 
  ggplot(aes(x = success)) +
  geom_density(fill = "hotpink", alpha = 0.5, color = "hotpink4", linewidth = 2) +
  xlab("posterior predictive # survived") +
  geom_vline(aes(xintercept = obs_success),
             color = "hotpink", linewidth = 2,
             linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(ppc_success, 0.025)), color = "hotpink4",
             linetype = "dotted") +
  geom_vline(aes(xintercept = quantile(ppc_success, 0.975)), color = "hotpink4",
             linetype = "dotted") -> ppc_surv_plot

# Predicted performance for fecundity model
predicted <- colMeans(y_rep_fecun)
observed <- cg_model_fecun$seed_count

give_me_R2 <- function(preds,actual){
  rss <- sum(( preds - actual ) ^ 2)  ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
  return(rsq)
}
give_me_R2(log(predicted), log(observed)) # 0.5812148

tibble(pred = log(predicted),
       obs = log(observed)) %>% 
  ggplot(aes(x = pred, y = obs)) +
  geom_point(shape = 1, stroke = 1.2, alpha = 0.1, size = 0.5) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", linewidth = 2, color = "gray47") +
  labs(x = "ln(predicted seed count)", y = "ln(observed seed count)") +
  annotate(geom = "text", x = 1.5, y = 10, label = 'R^2==0.581', parse=TRUE, size = 5) +
  lims(x = c(0,10), y = c(0,10))-> fecun_perform

# Predicted performance for survival model
confusion_data_vec <- function(pred_matrix, obs) {
  apply(pred_matrix, 1, function(pred) {
    n <- length(obs)
    true_pos  <- sum(pred == 1 & obs == 1) 
    true_neg  <- sum(pred == 0 & obs == 0) 
    false_pos <- sum(pred == 1 & obs == 0) 
    false_neg <- sum(pred == 0 & obs == 1) 
    c(true_pos = true_pos, true_neg = true_neg, 
      false_pos = false_pos, false_neg = false_neg)
  })
}

confusion_out <- confusion_data_vec(pred = y_rep_surv, obs = cg_model$survived)

# Compute summary statistics (mean and quantiles)
summary_conf <- apply(as.matrix(confusion_out), 1, function(x) {
  c(mean = mean(x),
    q2.5 = quantile(x, 0.025),
    q97.5 = quantile(x, 0.975))
})

# Convert to data frame for plotting
summary_df <- as.data.frame(t(summary_conf))
summary_df$Metric <- rownames(summary_df)

# Map confusion categories to 2x2 layout
plot_data <- data.frame(
  Metric = c("true_neg", "false_pos", "false_neg", "true_pos"),
  Row    = c(1, 1, 2, 2),  # y-axis position
  Col    = c(1, 2, 1, 2)   # x-axis position
)

# Merge with summary stats
plot_data <- merge(plot_data, summary_df, by = "Metric")

# Create label: mean and interval
plot_data$label <- sprintf(
  "%s\n%.0f",
  tolower(gsub("_", " ", plot_data$Metric)),
  plot_data$mean
)

# Ensure Row and Col are factors to treat grid as discrete space
plot_data$Col <- factor(plot_data$Col, levels = c(1, 2), labels = c("predicted: 0", "predicted: 1"))
plot_data$Row <- factor(plot_data$Row, levels = c(2, 1), labels = c("observed: 1", "observed: 0"))  # flip Y for top-down

plot_data$color <- c("bad", "bad", "good", "good")

ggplot(plot_data, aes(x = Col, y = Row)) +
  # Circles scaled by square root of mean
  geom_point(aes(size = mean, fill = color), shape = 21) +
  
  # Labels with mean and 95% interval
  geom_text(aes(label = label), size = 5) +
  
  coord_fixed() +
  scale_size_area(max_size = 60) +
  
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("gray67", "palegreen")) -> surv_perform

# Calculate accuracy 

accuracy <- function(tp, tn, n){
  (tp + tn) / n
}

accuracy(tp = plot_data %>% filter(Metric == "true_pos") %>% pull(mean),
         tn = plot_data %>% filter(Metric == "true_neg") %>% pull(mean),
         n = length(cg_model$survived)) 
# 0.663305

png("figs/FigS3_PPC.png", height = 7.3, width = 4.3, units = "in", res = 300)
ppc_mean_plot / ppc_sd_plot / ppc_surv_plot + plot_annotation(tag_level = "a")
dev.off()

png("figs/FigS4_modelperform.png", height = 4.5, width = 10, units = "in", res = 300)
surv_perform + fecun_perform + plot_annotation(tag_levels = "a")
dev.off()
