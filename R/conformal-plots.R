# conformal-plots.R
#
# Plots simulation results from conformal-runtime-sims.R.
# Produces boxplots and summary plots of coverage and interval length.
#-------------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggridges)
library(ggtext)
#-------------------------------------------------------------------------------
# Configuration: set input file here
sim_results_file <- '../output/sim_results.csv'
resdf <- read.csv(sim_results_file)
#-------------------------------------------------------------------------------
# Recoding lookups (used throughout)

method_labels <- c(
  dml       = "DML",
  wcp       = "Weighted",
  dml_naive = "DML, ignoring runtime confounding",
  wcp_naive = "WCP, ignoring runtime confounding",
  init      = "Init"
)

score_labels <- c(
  abs_resid = "Absolute residual score",
  quantile  = "Quantile score"
)

metric_labels <- c(
  cov = "Coverage Rate",
  len = "Interval Length"
)

#' Apply standard method/score/metric recoding to a long-format dataframe
recode_labels <- function(df) {
  df %>%
    mutate(
      method = recode(method, !!!method_labels),
      score  = recode(score, !!!score_labels),
      metric = recode(metric, !!!metric_labels)
    )
}
#-------------------------------------------------------------------------------
# Reshape to long format
# Column naming from sims: {method}_{metric} where metric is cov or len

resdf_long <- resdf %>%
  pivot_longer(
    cols = matches("_(cov|len)$"),
    names_to = c("method", "metric"),
    names_pattern = "(dml_naive|wcp_naive|dml|wcp|init)_(cov|len)",
    values_to = "value"
  )
#-------------------------------------------------------------------------------
# Main plot

resdf_long %>%
  filter(method != 'init', d2 == 10) %>%
  recode_labels() %>%
  mutate(samp = as.factor(n)) %>%
  ggplot() +
  geom_boxplot(aes(x = samp, y = value, fill = method), alpha = 0.5) +
  facet_grid(metric ~ score, scales = 'free_y') +
  geom_hline(data = . %>% filter(metric == 'Coverage Rate') %>% mutate(ref = 0.9),
             aes(yintercept = ref),
             linetype = "dashed", color = "red", linewidth = 1) +
  theme_bw() +
  theme(legend.position = 'bottom',
        title = element_markdown(),
        legend.text = element_markdown()) +
  labs(title = '**Conformal prediction of counterfactuals in target population**',
       subtitle = 'Average coverage and length of prediction intervals, 10 runtime confounders',
       x = 'Overall sample size, n',
       y = 'Operating Characteristic',
       fill = 'Method')
ggsave('../figures/coverage_length_conformal_target-main.pdf', width = 8.3, height = 5, units = 'in')
#-------------------------------------------------------------------------------
# Main results, by A

resdf_long %>%
  filter(method != 'init', d2 == 10) %>%
  recode_labels() %>%
  mutate(samp = as.factor(n), A_lab = paste0('A=', A)) %>%
  ggplot() +
  geom_boxplot(aes(x = samp, y = value, fill = method), alpha = 0.5) +
  facet_grid(metric + A_lab ~ score, scales = 'free_y') +
  geom_hline(data = . %>% filter(metric == 'Coverage Rate') %>% mutate(ref = 0.9),
             aes(yintercept = ref),
             linetype = "dashed", color = "red", linewidth = 1) +
  theme_bw() +
  theme(legend.position = 'bottom',
        title = element_markdown(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown()) +
  labs(title = '**Conformal prediction of counterfactuals in target population**',
       subtitle = 'Average coverage and length of prediction intervals',
       x = 'Overall sample size, n',
       y = 'Operating Characteristic',
       fill = 'Method')
ggsave('../figures/coverage_length_conformal_target-byA.pdf', width = 8, height = 8, units = 'in')
#-------------------------------------------------------------------------------
# Summary of length and coverage by n and method

plotdf <- resdf_long %>%
  filter(method != 'init', rho_val == 0.9) %>%
  recode_labels() %>%
  group_by(n, method, metric, d2, score) %>%
  summarise(
    avg_value = mean(value),
    num       = n(),
    sd_value  = sd(value),
    .groups   = 'drop'
  )

plotdf %>% filter(n == 5000) %>%
  ggplot(aes(x = d2, y = avg_value, color = method)) +
  geom_point(size = 1.8, fill = 'white', shape = 6, stroke = 1.2, alpha = 0.8,
             position = position_dodge(1)) +
  geom_errorbar(aes(ymin = avg_value - sd_value,
                    ymax = avg_value + sd_value),
                width = 0.2, position = position_dodge(1)) +
  geom_hline(data = subset(plotdf, metric == 'Coverage Rate'),
             aes(yintercept = 0.9),
             linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(breaks = c(5, 10, 15)) +
  facet_grid(metric ~ score, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'bottom',
        title = element_markdown(),
        legend.text = element_markdown()) +
  labs(title = '**Conformal prediction of counterfactuals in target population**',
       subtitle = 'Average coverage and length, fixing n=5000 and varying number of runtime confounders',
       x = 'Total number of runtime confounders',
       y = 'Operating Characteristic',
       color = 'Method')
ggsave('../figures/coverage_over_d2.pdf', width = 8.3, height = 5, units = 'in')
#-------------------------------------------------------------------------------
# Appendix plots: vary runtime confounding dimension

for (dd in c(5, 15)) {

  resdf_long %>%
    filter(method != 'init', d2 == dd, rho_val == 0.9) %>%
    recode_labels() %>%
    mutate(samp = as.factor(n), A_lab = paste0('A=', A)) %>%
    ggplot() +
    geom_boxplot(aes(x = samp, y = value, fill = method), alpha = 0.5) +
    facet_grid(metric + A_lab ~ score, scales = 'free_y') +
    geom_hline(data = . %>% filter(metric == 'Coverage Rate') %>% mutate(ref = 0.9),
               aes(yintercept = ref),
               linetype = "dashed", color = "red") +
    theme_bw() +
    theme(legend.position = 'bottom',
          title = element_markdown(),
          legend.text = element_markdown()) +
    labs(title = '**Conformal prediction of counterfactuals in target population**',
         subtitle = paste0('Average coverage and length of prediction intervals, ', dd, ' runtime confounders'),
         x = 'Overall sample size, n',
         y = 'Operating Characteristic',
         fill = 'Method')
  ggsave(paste0('../figures/coverage_length_conformal_target-main-', dd, '.pdf'),
         width = 8.3, height = 7, units = 'in')
}

# Summary across all d2 values
plotdf %>%
  ggplot(aes(x = n, y = avg_value, color = method)) +
  geom_point(size = 1.8, fill = 'white', shape = 6, stroke = 1.2, alpha = 0.8,
             position = position_dodge(500)) +
  geom_errorbar(aes(ymin = avg_value - sd_value,
                    ymax = avg_value + sd_value),
                width = 0.2, position = position_dodge(500)) +
  geom_hline(data = subset(plotdf, metric == 'Coverage Rate'),
             aes(yintercept = 0.9),
             linetype = "dashed", color = "red") +
  facet_grid(metric + factor(d2) ~ score, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'bottom',
        title = element_markdown(),
        legend.text = element_markdown()) +
  labs(title = '**Conformal prediction of counterfactuals in target population**',
       subtitle = 'Average coverage and length, varying n and number of runtime confounders',
       x = 'Overall sample size, n',
       y = 'Operating Characteristic',
       color = 'Method')
#-------------------------------------------------------------------------------
# Appendix: varying share of source data

resdf_long %>%
  filter(d2 == 10, method != 'init', n == 5000) %>%
  recode_labels() %>%
  mutate(A_lab = paste0('A=', A)) %>%
  ggplot() +
  geom_boxplot(aes(x = factor(rho_val), y = value, fill = method), alpha = 0.5) +
  facet_grid(metric + A_lab ~ score, scales = 'free_y') +
  geom_hline(data = . %>% filter(metric == 'Coverage Rate') %>% mutate(ref = 0.9),
             aes(yintercept = ref),
             linetype = "dashed", color = "red", linewidth = 1) +
  theme_bw() +
  theme(legend.position = 'bottom',
        title = element_markdown(),
        legend.text = element_markdown()) +
  labs(title = '**Conformal prediction of counterfactuals in target population**',
       subtitle = 'Varying share of source data, fixing 10 runtime confounders and n=5000',
       x = 'Share of source data, P(S=1)',
       y = 'Operating Characteristic',
       fill = 'Method')
ggsave('../figures/coverage_length_conformal_target-rhovary.pdf', width = 10.3, height = 7, units = 'in')
