# plotting-acic.R
#
# Plots results from the semi-synthetic ACIC data application (clean-acic.R).
#-------------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggtext)
#-------------------------------------------------------------------------------
# Configuration: set input file here
acic_results_file <- '../output/acic-results.csv'
res <- read.csv(acic_results_file)
#-------------------------------------------------------------------------------
# Recoding lookups

method_labels <- c(
  dml       = "DML",
  dml_naive = "DML, ignoring runtime confounding",
  wcp       = "Weighted",
  wcp_naive = "Weighted (ignoring confounding)"
)

score_labels <- c(
  abs_resid = "Absolute Residual",
  quantile  = "Quantile"
)

metric_labels <- c(
  cov = "Coverage Rate",
  len = "Interval Length"
)
#-------------------------------------------------------------------------------
# Reshape to long format
# Column naming from clean-acic.R: {method}_{metric}

reslong <- res %>%
  pivot_longer(
    cols = matches("_(cov|len)$"),
    names_to = c("method", "metric"),
    names_pattern = "(dml_naive|wcp_naive|dml|wcp)_(cov|len)",
    values_to = "value"
  ) %>%
  filter(method != 'wcp_naive') %>%
  mutate(
    method = recode(method, !!!method_labels),
    metric = recode(metric, !!!metric_labels),
    score  = recode(score, !!!score_labels),
    A_lab  = paste0('A=', A)
  )
#-------------------------------------------------------------------------------
#' Helper to create a standard ACIC boxplot
#'
#' @param data      Filtered long-format dataframe
#' @param subtitle  Plot subtitle
#' @param facet_formula Facet formula (e.g., metric ~ score)
make_acic_plot <- function(data, subtitle, facet_formula = metric ~ score) {
  data %>%
    ggplot(aes(x = factor(nn), y = value, fill = method)) +
    geom_boxplot(position = 'dodge', alpha = 0.5) +
    facet_grid(facet_formula, scales = 'free') +
    geom_hline(data = . %>% filter(metric == 'Coverage Rate') %>% mutate(ref = 0.9),
               aes(yintercept = ref),
               linetype = "dashed", color = "red", linewidth = 1) +
    theme_bw() +
    theme(legend.position = 'bottom',
          title = element_markdown(),
          legend.text = element_markdown()) +
    labs(
      title = "**Conformal prediction of counterfactuals in target population**",
      subtitle = subtitle,
      x = "Overall sample size, n",
      y = "Operating Characteristic",
      fill = 'Method'
    )
}
#-------------------------------------------------------------------------------
# Baseline plot (moderate runtime confounding)

make_acic_plot(
  reslong %>% filter(rc_sev == 'moderate'),
  subtitle = 'Semi-synthetic ACIC data'
) +
  geom_blank(data = reslong %>% filter(metric == "Coverage Rate"), aes(y = 0.65)) +
  geom_blank(data = reslong %>% filter(metric == "Coverage Rate"), aes(y = 1))
ggsave('../figures/acic-baseline.pdf', width = 8.3, height = 5)
#-------------------------------------------------------------------------------
# Baseline, separated by A

make_acic_plot(
  reslong %>% filter(rc_sev == 'moderate'),
  subtitle = 'Semi-synthetic ACIC data, results stratified by counterfactual outcome a=0,1',
  facet_formula = metric + A_lab ~ score
)
ggsave('../figures/acic-baseline-byA.pdf', width = 8, height = 8)
#-------------------------------------------------------------------------------
# Mild runtime confounding

make_acic_plot(
  reslong %>% filter(rc_sev == 'mild'),
  subtitle = 'Semi-synthetic ACIC data, mild runtime confounding'
)
ggsave('../figures/acic-mild.pdf', width = 8, height = 6)
#-------------------------------------------------------------------------------
# Severe runtime confounding

make_acic_plot(
  reslong %>% filter(rc_sev == 'severe'),
  subtitle = 'Semi-synthetic ACIC data, severe runtime confounding'
)
ggsave('../figures/acic-severe.pdf', width = 8, height = 6)
