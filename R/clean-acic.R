# clean-acic.R
#
# Semi-synthetic data application using the 2018 ACIC competition data,
# following the DGP outlined in Lei and Candes (2021). Compares:
#   - DML conformal (with and without runtime confounding adjustment)
#   - Weighted conformal (with and without runtime confounding adjustment)
# across varying sample sizes, score functions, and severity of runtime
# confounding.
#
# Output: CSV in ../output/ with per-replicate coverage and interval lengths.
#-------------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(foreach)
library(doParallel)
library(randomForest)
library(SuperLearner)
source('estimation.R')
#-------------------------------------------------------------------------------
# Load in data + emulate setting in Lei and Candes (2021)

acic_data <- read.csv('../data/synthetic_data.csv') %>%
  mutate(S3 = factor(S3),
         C1 = factor(C1),
         C2 = factor(C2),
         C3 = factor(C3),
         XC = factor(XC),
         schoolid = factor(schoolid)) %>%
  mutate(carvalho_tau = 0.228 + 0.05 * as.numeric(X1 < 0.07) -
           0.05 * as.numeric(X2 < -0.69) -
           0.08 * as.numeric(C1 %in% c(1, 13, 14)))

# Split into training and generation sets
idx <- sample(1:nrow(acic_data), 2079, replace = FALSE)
acic_data$train <- as.numeric(1:nrow(acic_data) %in% idx)

# Set up model matrices
all_covs <- c('S3', 'C1', 'C2', 'C3', 'XC', 'X1', 'X2', 'X3', 'X4', 'X5')

covs    <- acic_data %>% filter(train == 1) %>% select(all_of(all_covs))
covsZ   <- acic_data %>% filter(train == 1) %>% select(all_of(c(all_covs, 'Z')))
covsY   <- acic_data %>% filter(train == 1) %>% select(all_of(c(all_covs, 'Y')))
covs0   <- acic_data %>% filter(Z == 0, train == 1) %>% select(all_of(all_covs))
covs1   <- acic_data %>% filter(Z == 1, train == 1) %>% select(all_of(all_covs))
covs_gen <- acic_data %>% filter(train == 0) %>% select(all_of(all_covs))

# Fit model for E[Y(0)]
m0mod <- randomForest(Y ~ ., data = covsY, ntree = 500, importance = TRUE)

# Fit P(A=1|X)
gmod <- randomForest(factor(Z) ~ ., data = covsZ, ntree = 500, importance = TRUE)

# Fit 25th and 75th percentiles of observed Y0 and Y1
Y0 <- acic_data %>% filter(Z == 0, train == 1) %>% pull(Y)
Y1 <- acic_data %>% filter(Z == 1, train == 1) %>% pull(Y)
x0 <- model.matrix(~ . - 1, data = covs0)
x1 <- model.matrix(~ . - 1, data = covs1)

r0mod <- grf::quantile_forest(x0, Y0, quantiles = c(0.25, 0.75))
r1mod <- grf::quantile_forest(x1, Y1, quantiles = c(0.25, 0.75))

carvalho_tau <- function(covs) {
  0.228 + 0.05 * as.numeric(covs$X1 < 0.07) -
    0.05 * as.numeric(covs$X2 < -0.69) -
    0.08 * as.numeric(covs$C1 %in% c(1, 13, 14))
}

# All covariates and partitions into always-observed (V) and runtime (W)
all_vars <- c('X1', 'X2', 'X3', 'X4', 'X5', 'S3', 'C1', 'C2', 'C3', 'XC')
#-------------------------------------------------------------------------------
#' Simulate semi-synthetic data from ACIC generating process
#'
#' @param covs       Covariate dataframe to resample from
#' @param mu0mod     Random forest for E[Y(0)|X]
#' @param gmod       Random forest for P(A=1|X)
#' @param r0mod      Quantile forest for IQR of Y(0)
#' @param r1mod      Quantile forest for IQR of Y(1)
#' @param Smod       Random forest for source assignment model
#' @param bval       Intercept for source assignment (passed explicitly)
#' @param nsim       Number of observations to generate
sim_data <- function(covs, mu0mod, gmod, r0mod, r1mod, Smod,
                     bval, nsim = 3000) {

  idx <- sample(1:nrow(covs), nsim, replace = TRUE)
  X <- covs[idx, ]
  Xmat <- model.matrix(~ . - 1, data = X)

  # Generate S
  S_lin <- qlogis(pmin(pmax(predict(Smod, newdata = X, type = 'prob')[, 2], 0.05), 0.95))
  Sprobs <- plogis(bval - S_lin)
  S <- rbinom(nrow(X), 1, Sprobs)

  # Generate treatment
  Aprobs <- predict(gmod, newdata = X, type = 'prob')[, 2]
  Aprobs <- pmax(pmin(Aprobs, 0.9), 0.1)
  A <- rbinom(nsim, 1, Aprobs)

  # Generate Y0
  IQR0 <- predict(r0mod, newdata = Xmat, quantiles = c(0.25, 0.75))$predictions
  IQR0 <- IQR0[, 2] - IQR0[, 1]
  mhat <- predict(mu0mod, newdata = X)
  Y0 <- mhat + 0.5 * IQR0 * rnorm(nsim)

  # Generate Y1
  IQR1 <- predict(r1mod, newdata = Xmat, quantiles = c(0.25, 0.75))$predictions
  IQR1 <- IQR1[, 2] - IQR1[, 1]
  Y1 <- mhat + carvalho_tau(X) + 0.5 * IQR1 * rnorm(nsim)

  # Observed Y
  Y <- A * Y1 + (1 - A) * Y0

  cbind(data.frame(Y = Y, A = A, S = S, Y1 = Y1, Y0 = Y0), X)
}
#-------------------------------------------------------------------------------
#' Compute coverage and interval length for target population
eval_coverage <- function(res_df, truth) {
  target <- res_df %>% filter(S == 0)
  truth_target <- truth[res_df$S == 0]
  list(cov = mean((truth_target >= target$int_L) & (truth_target <= target$int_U)),
       len = mean(target$int_U - target$int_L))
}
#-------------------------------------------------------------------------------
#' Run one ACIC simulation replicate
#'
#' @param the_data   Simulated data from sim_data()
#' @param vars_v     Names of always-observed covariates
#' @param vars_w     Names of runtime confounders
#' @param score      Conformal score type
#' @param alpha      Nominal miscoverage level
#' @param la         Named list of learner libraries
run_one_acic_sim <- function(the_data, vars_v, vars_w, score, alpha, la) {

  Y    <- the_data$Y
  A    <- the_data$A
  S    <- the_data$S
  V_df <- data.frame(model.matrix(~ 0 + ., data = the_data[, vars_v, drop = FALSE]))
  W_df <- data.frame(model.matrix(~ 0 + ., data = the_data[, vars_w, drop = FALSE]))

  # --- Weighted conformal, adjusting for runtime confounding ---
  wcp_a1 <- wcp_rc(Y, A, V_df, W_df, S,
                    score = score, mu_learners = la$mu,
                    g_learners = la$g, kappa_learners = la$kappa,
                    alpha = alpha)
  wcp_a1$Y1 <- the_data$Y1[S == 0]
  wcp_a1_eval <- eval_coverage(wcp_a1, the_data$Y1)

  wcp_a0 <- wcp_rc(Y, 1 - A, V_df, W_df, S,
                    score = score, mu_learners = la$mu,
                    g_learners = la$g, kappa_learners = la$kappa,
                    alpha = alpha)
  wcp_a0$Y0 <- the_data$Y0[S == 0]
  wcp_a0_eval <- eval_coverage(wcp_a0, the_data$Y0)

  # --- Weighted conformal, ignoring runtime confounding ---
  wcp_naive_a1 <- wcp_rc(Y, A, V_df, NULL, S,
                          score = score, mu_learners = la$mu,
                          g_learners = la$g, kappa_learners = la$kappa,
                          alpha = alpha)
  wcp_naive_a1$Y1 <- the_data$Y1[S == 0]
  wcp_naive_a1_eval <- eval_coverage(wcp_naive_a1, the_data$Y1)

  wcp_naive_a0 <- wcp_rc(Y, 1 - A, V_df, NULL, S,
                          score = score, mu_learners = la$mu,
                          g_learners = la$g, kappa_learners = la$kappa,
                          alpha = alpha)
  wcp_naive_a0$Y0 <- the_data$Y0[S == 0]
  wcp_naive_a0_eval <- eval_coverage(wcp_naive_a0, the_data$Y0)

  # --- DML conformal, adjusting for runtime confounding ---
  dml_a1 <- est_r_dml(Y, A, V_df, W_df, S,
                       score = score, mu_learners = la$mu,
                       g_learners = la$g, kappa_learners = la$kappa,
                       m_learners = la$m, q_learners = la$q,
                       alpha = alpha) %>%
    mutate(Y1 = the_data$Y1, Y0 = the_data$Y0)
  dml_a1_eval <- eval_coverage(dml_a1, the_data$Y1)

  dml_a0 <- est_r_dml(Y, 1 - A, V_df, W_df, S,
                       score = score, mu_learners = la$mu,
                       g_learners = la$g, kappa_learners = la$kappa,
                       m_learners = la$m, q_learners = la$q,
                       alpha = alpha) %>%
    mutate(Y1 = the_data$Y1, Y0 = the_data$Y0)
  dml_a0_eval <- eval_coverage(dml_a0, the_data$Y0)

  # --- DML conformal, ignoring runtime confounding ---
  dml_naive_a1 <- est_r_dml(Y, A, V_df, NULL, S,
                             score = score, mu_learners = la$mu,
                             g_learners = la$g, kappa_learners = la$kappa,
                             m_learners = la$m, q_learners = la$q,
                             alpha = alpha) %>%
    mutate(Y1 = the_data$Y1, Y0 = the_data$Y0)
  dml_naive_a1_eval <- eval_coverage(dml_naive_a1, the_data$Y1)

  dml_naive_a0 <- est_r_dml(Y, 1 - A, V_df, NULL, S,
                             score = score, mu_learners = la$mu,
                             g_learners = la$g, kappa_learners = la$kappa,
                             m_learners = la$m, q_learners = la$q,
                             alpha = alpha) %>%
    mutate(Y1 = the_data$Y1, Y0 = the_data$Y0)
  dml_naive_a0_eval <- eval_coverage(dml_naive_a0, the_data$Y0)

  # --- Combine results ---
  # Column naming: {method}_{metric}, consistent with conformal-runtime-sims.R
  df_a1 <- data.frame(
    dml_cov = dml_a1_eval$cov,             dml_naive_cov = dml_naive_a1_eval$cov,
    wcp_cov = wcp_a1_eval$cov,             wcp_naive_cov = wcp_naive_a1_eval$cov,
    dml_len = dml_a1_eval$len,             dml_naive_len = dml_naive_a1_eval$len,
    wcp_len = wcp_a1_eval$len,             wcp_naive_len = wcp_naive_a1_eval$len,
    A = 1)
  df_a0 <- data.frame(
    dml_cov = dml_a0_eval$cov,             dml_naive_cov = dml_naive_a0_eval$cov,
    wcp_cov = wcp_a0_eval$cov,             wcp_naive_cov = wcp_naive_a0_eval$cov,
    dml_len = dml_a0_eval$len,             dml_naive_len = dml_naive_a0_eval$len,
    wcp_len = wcp_a0_eval$len,             wcp_naive_len = wcp_naive_a0_eval$len,
    A = 0)

  rbind(df_a1, df_a0)
}
#-------------------------------------------------------------------------------
# Set up params for semi-synthetic data

nsim   <- 500
maxits <- 3

rf_bin <- create.Learner("SL.ranger",
                         params = list(probability      = TRUE,
                                       num.trees        = 350,
                                       max.depth        = 0,
                                       min.node.size    = 5,
                                       replace          = FALSE,
                                       sample.fraction  = 0.7,
                                       splitrule        = "hellinger"),
                         name_prefix = "ranger_bin")

rf_cont <- create.Learner("SL.ranger",
                          params = list(num.trees        = 350,
                                        max.depth        = 0,
                                        min.node.size    = 5,
                                        replace          = FALSE,
                                        sample.fraction  = 0.7,
                                        splitrule        = "variance"),
                          name_prefix = "ranger_cont")

g_learners <- kappa_learners <- q_learners <- c('SL.glm', 'SL.glmnet', 'ranger_bin_1')
mu_learners <- m_learners <- c('SL.glm', 'SL.glmnet', 'ranger_cont_1')
alpha <- 0.1

learner_args <- list(mu = mu_learners, g = g_learners,
                     kappa = kappa_learners, m = m_learners, q = q_learners)

param_grid <- expand.grid(
  score  = c("abs_resid", "quantile"),
  nn     = c(2500, 5000, 7500, 1e4),
  rc_sev = c('mild', 'severe'),
  stringsAsFactors = FALSE
)
#-------------------------------------------------------------------------------
# Main simulation

doParallel::registerDoParallel(cores = 40)

resdf <- data.frame()
for (pp in 1:nrow(param_grid)) {

  print(pp)
  score  <- param_grid$score[pp]
  nn     <- param_grid$nn[pp]
  rc_sev <- param_grid$rc_sev[pp]

  # Determine V/W partition based on severity of runtime confounding
  if (rc_sev == 'moderate') {
    vars_v <- c('X1', 'X2', 'X3', 'X4', 'X5', 'XC')
  } else if (rc_sev == 'severe') {
    vars_v <- c('X3', 'X4', 'X5', 'XC')
  } else if (rc_sev == 'mild') {
    vars_v <- c('X1', 'X2', 'X3', 'X4', 'X5', 'XC', 'S3', 'C1')
  }
  vars_w <- setdiff(all_vars, vars_v)

  # Create generating model for S, using transformation of treatment propensity
  Smod <- randomForest(factor(Z) ~ .,
                       data = covsZ %>% select(all_of(c('Z', vars_v))),
                       ntree = 500)

  # Find intercept b so that P(S=1) = 0.9
  cal_idx <- sample(1:nrow(covs_gen), 1e5, replace = TRUE)
  X_cal   <- covs_gen[cal_idx, ]
  S_lin   <- qlogis(pmin(pmax(predict(Smod, newdata = X_cal, type = 'prob')[, 2], 0.05), 0.95))
  s_root  <- function(b) mean(plogis(-S_lin + b)) - 0.9
  bval    <- uniroot(s_root, c(-10, 10))$root

  curr_res <- foreach(ss = 1:nsim, .combine = rbind) %dopar% {
    out <- NULL
    for (it in 1:maxits) {
      attempt <- try({
        the_data <- sim_data(covs_gen, m0mod, gmod, r0mod, r1mod, Smod,
                             bval = bval, nsim = nn)
        run_one_acic_sim(the_data, vars_v, vars_w, score, alpha, learner_args)
      }, silent = TRUE)
      if (!inherits(attempt, "try-error")) { out <- attempt; break }
    }
    if (is.null(out)) stop("All retries failed for this iteration")
    out
  }

  # Store params for this grid point
  curr_res$score  <- param_grid$score[pp]
  curr_res$nn     <- param_grid$nn[pp]
  curr_res$rc_sev <- param_grid$rc_sev[pp]

  resdf <- rbind(resdf, curr_res)
}

doParallel::stopImplicitCluster()
#-------------------------------------------------------------------------------
# Save results


write.csv(resdf,
          '../output/acic-results.csv',
          row.names = FALSE)
