# conformal-runtime-sims.R
#
# Simulation study for debiased conformal prediction of counterfactual
# outcomes under runtime confounding. Compares:
#   - DML conformal (proposed method, with and without runtime confounding adj.)
#   - Weighted conformal (with and without runtime confounding adj.)
# across varying sample sizes, score functions, and numbers of runtime confounders.
#
# Output: CSV in ../output/ with per-replicate coverage and interval lengths.
#-------------------------------------------------------------------------------
rm(list=ls())

library(tidyverse)
library(SuperLearner)
library(foreach)
library(doParallel)
source('estimation.R')
#-------------------------------------------------------------------------------
# Fixed simulation parameters

alpha <- 0.1                # nominal miscoverage level
p1 <- 15                    # dimension of always-observed covariates V
d1 <- 5                     # number of active dimensions in V
p2 <- 15                    # dimension of runtime confounders W

# SuperLearner libraries
g_learners <- mu_learners <- kappa_learners <- m_learners <- q_learners <-
  c('SL.glm', 'SL.glmnet', 'SL.ranger')
#-------------------------------------------------------------------------------
# Parameters that vary across simulation settings

n_grid     <- 5000                          # sample sizes
score_grid <- c('abs_resid', 'quantile')    # conformal score functions
dw_grid    <- c(5, 15)                      # number of active dims in W
rho_grid   <- c(0.9)                        # P(S=1), share of source data

param_combos <- expand.grid(n = n_grid, score = score_grid,
                            d2 = dw_grid, rho_grid = rho_grid)
#-------------------------------------------------------------------------------
#' Find intercept b so that P(S=1) = targ under the source assignment model
get_S_b <- function(V, d1, targ = 0.8) {
  init_comp <- -1/sqrt(d1) * rowSums(V[, 1:d1])
  prob_func <- function(b) mean(plogis(b + init_comp)) - targ
  uniroot(prob_func, c(-50, 50))$root
}
#-------------------------------------------------------------------------------
#' Compute coverage and interval length for a set of conformal results
#'
#' @param res_df Results dataframe from wcp_rc or est_r_dml (filtered to S==0)
#' @param truth  Vector of true potential outcomes for the target population
#' @return Named list with coverage and length
eval_coverage <- function(res_df, truth) {
  covered <- (truth >= res_df$int_L) & (truth <= res_df$int_U)
  list(cov = mean(covered),
       len = mean(res_df$int_U - res_df$int_L))
}
#-------------------------------------------------------------------------------
#' Run a single simulation replicate
#'
#' @param n         Sample size
#' @param score     Conformal score type ('abs_resid' or 'quantile')
#' @param d2        Number of active dimensions in W
#' @param b_shift   Intercept for source assignment model
#' @param p1,p2,d1  Dimension parameters
#' @param alpha     Nominal miscoverage level
#' @param learner_args Named list of learner libraries
#' @return One-row (two rows: A=1, A=0) dataframe with coverage/length results
run_one_sim <- function(n, score, d2, b_shift, p1, p2, d1, alpha,
                        learner_args) {

  # --- Simulate data ---
  V <- mvtnorm::rmvnorm(n, mean = rep(0, p1), sigma = diag(p1))
  W <- mvtnorm::rmvnorm(n, mean = rep(0, p2), sigma = diag(p2))

  # Source vs target and treatment assignment
  S <- rbinom(n, 1, plogis(b_shift - 1/sqrt(d1) * rowSums(V[, 1:d1])))
  A <- rbinom(n, 1, plogis(1/sqrt(d1 + d2) *
                              (rowSums(V[, 1:d1]) - 2 * rowSums(W[, 1:d2]))))

  # Potential outcomes
  mu_Y <- d1/(d1 + d2) * (rowSums(V[, 1:d1]) + 2 * rowSums(W[, 1:d2]))
  sd_Y <- d1/(d1 + d2) * (rowSums(V[, 1:d1]) + rowSums(W[, 1:d2]))
  Y0 <- mu_Y + rnorm(n, sd = abs(sd_Y))
  Y1 <- mu_Y + rnorm(n, sd = abs(sd_Y))
  Y  <- ifelse(A == 0, Y0, Y1)

  # Covariate dataframes
  V_df <- data.frame(V); colnames(V_df) <- paste0("V", seq_len(ncol(V_df)))
  W_df <- data.frame(W); colnames(W_df) <- paste0("W", seq_len(ncol(W_df)))

  la <- learner_args  # shorthand

  # --- Weighted conformal, adjusting for runtime confounding ---
  wcp_a1 <- wcp_rc(Y, A, V_df, W_df, S,
                    score = score, mu_learners = la$mu,
                    g_learners = la$g, kappa_learners = la$kappa,
                    alpha = alpha)
  wcp_a1$Y1 <- Y1[S == 0]
  wcp_a0 <- wcp_rc(Y, 1 - A, V_df, W_df, S,
                    score = score, mu_learners = la$mu,
                    g_learners = la$g, kappa_learners = la$kappa,
                    alpha = alpha)
  wcp_a0$Y0 <- Y0[S == 0]

  wcp_a1_eval <- eval_coverage(wcp_a1 %>% filter(S == 0), wcp_a1$Y1[wcp_a1$S == 0])
  wcp_a0_eval <- eval_coverage(wcp_a0 %>% filter(S == 0), wcp_a0$Y0[wcp_a0$S == 0])

  # --- DML conformal, adjusting for runtime confounding ---
  dml_a1 <- est_r_dml(Y, A, V_df, W_df, S,
                       score = score, mu_learners = la$mu,
                       g_learners = la$g, kappa_learners = la$kappa,
                       m_learners = la$m, q_learners = la$q,
                       alpha = alpha) %>%
    mutate(Y0 = Y0, Y1 = Y1)
  dml_a0 <- est_r_dml(Y, 1 - A, V_df, W_df, S,
                       score = score, mu_learners = la$mu,
                       g_learners = la$g, kappa_learners = la$kappa,
                       m_learners = la$m, q_learners = la$q,
                       alpha = alpha) %>%
    mutate(Y0 = Y0, Y1 = Y1)

  dml_a1_main <- dml_a1 %>% filter(S == 0)
  dml_a0_main <- dml_a0 %>% filter(S == 0)

  dml_a1_eval  <- eval_coverage(dml_a1_main, dml_a1_main$Y1)
  dml_a0_eval  <- eval_coverage(dml_a0_main, dml_a0_main$Y0)

  # Initial (pre-debiasing) DML coverage
  init_a1_eval <- list(
    cov = mean((dml_a1_main$Y1 >= dml_a1_main$int_L_init) &
               (dml_a1_main$Y1 <= dml_a1_main$int_U_init)),
    len = mean(dml_a1_main$int_U_init - dml_a1_main$int_L_init))
  init_a0_eval <- list(
    cov = mean((dml_a0_main$Y0 >= dml_a0_main$int_L_init) &
               (dml_a0_main$Y0 <= dml_a0_main$int_U_init)),
    len = mean(dml_a0_main$int_U_init - dml_a0_main$int_L_init))

  # --- DML ignoring runtime confounding (W=NULL) ---
  dml_naive_a1 <- est_r_dml(Y, A, V_df, NULL, S,
                             score = score, mu_learners = la$mu,
                             g_learners = la$g, kappa_learners = la$kappa,
                             m_learners = la$m, q_learners = la$q,
                             alpha = alpha) %>%
    mutate(Y0 = Y0, Y1 = Y1)
  dml_naive_a0 <- est_r_dml(Y, 1 - A, V_df, NULL, S,
                             score = score, mu_learners = la$mu,
                             g_learners = la$g, kappa_learners = la$kappa,
                             m_learners = la$m, q_learners = la$q,
                             alpha = alpha) %>%
    mutate(Y0 = Y0, Y1 = Y1)

  dml_naive_a1_eval <- eval_coverage(dml_naive_a1 %>% filter(S == 0),
                                     dml_naive_a1$Y1[dml_naive_a1$S == 0])
  dml_naive_a0_eval <- eval_coverage(dml_naive_a0 %>% filter(S == 0),
                                     dml_naive_a0$Y0[dml_naive_a0$S == 0])

  # --- Weighted conformal, ignoring runtime confounding (W=NULL) ---
  wcp_naive_a1 <- wcp_rc(Y, A, V_df, NULL, S,
                          score = score, mu_learners = la$mu,
                          g_learners = la$g, kappa_learners = la$kappa,
                          alpha = alpha)
  wcp_naive_a1$Y1 <- Y1[S == 0]
  wcp_naive_a0 <- wcp_rc(Y, 1 - A, V_df, NULL, S,
                          score = score, mu_learners = la$mu,
                          g_learners = la$g, kappa_learners = la$kappa,
                          alpha = alpha)
  wcp_naive_a0$Y0 <- Y0[S == 0]

  wcp_naive_a1_eval <- eval_coverage(wcp_naive_a1 %>% filter(S == 0),
                                     wcp_naive_a1$Y1[wcp_naive_a1$S == 0])
  wcp_naive_a0_eval <- eval_coverage(wcp_naive_a0 %>% filter(S == 0),
                                     wcp_naive_a0$Y0[wcp_naive_a0$S == 0])

  # --- Combine results ---
  # Column naming convention: {method}_{metric}
  #   method: dml, wcp, dml_naive, wcp_naive, init
  #   metric: cov (coverage), len (interval length)
  df_a1 <- data.frame(
    dml_cov = dml_a1_eval$cov,       dml_len = dml_a1_eval$len,
    wcp_cov = wcp_a1_eval$cov,       wcp_len = wcp_a1_eval$len,
    dml_naive_cov = dml_naive_a1_eval$cov,  dml_naive_len = dml_naive_a1_eval$len,
    wcp_naive_cov = wcp_naive_a1_eval$cov,  wcp_naive_len = wcp_naive_a1_eval$len,
    init_cov = init_a1_eval$cov,     init_len = init_a1_eval$len,
    n = n, score = score, A = 1, d2 = d2, alpha = alpha, rho_val = rho_val)
  df_a0 <- data.frame(
    dml_cov = dml_a0_eval$cov,       dml_len = dml_a0_eval$len,
    wcp_cov = wcp_a0_eval$cov,       wcp_len = wcp_a0_eval$len,
    dml_naive_cov = dml_naive_a0_eval$cov,  dml_naive_len = dml_naive_a0_eval$len,
    wcp_naive_cov = wcp_naive_a0_eval$cov,  wcp_naive_len = wcp_naive_a0_eval$len,
    init_cov = init_a0_eval$cov,     init_len = init_a0_eval$len,
    n = n, score = score, A = 0, d2 = d2, alpha = alpha, rho_val = rho_val)

  rbind(df_a1, df_a0)
}
#-------------------------------------------------------------------------------
# Main simulation loop

nsim <- 500
maxits <- 3     # max retries per replicate (SL occasionally fails numerically)

registerDoParallel(cores = 50)

learner_args <- list(mu = mu_learners, g = g_learners,
                     kappa = kappa_learners, m = m_learners, q = q_learners)

resdf <- data.frame()

for (ii in 1:nrow(param_combos)) {
  print(ii)

  n       <- param_combos[ii, 'n']
  score   <- param_combos[ii, 'score']
  d2      <- param_combos[ii, 'd2']
  rho_val <- param_combos[ii, 'rho_grid']

  # Calibrate intercept so P(S=1) = rho_val
  V_cal   <- mvtnorm::rmvnorm(1e6, mean = rep(0, p1), sigma = diag(p1))
  b_shift <- get_S_b(V_cal, d1, targ = rho_val)

  curr_res <- foreach(ss = 1:nsim, .combine = rbind,
                      .packages = c("dplyr")) %dopar% {
    out <- NULL
    for (it in 1:maxits) {
      attempt <- try({
        run_one_sim(n, score, d2, b_shift, p1, p2, d1, alpha, learner_args)
      }, silent = TRUE)
      if (!inherits(attempt, "try-error")) { out <- attempt; break }
    }
    if (is.null(out)) stop("All retries failed for this iteration")
    out
  }

  resdf <- rbind(resdf, curr_res)
}
#-------------------------------------------------------------------------------
# Write results

resdf <- resdf %>%
  mutate(
    mu_learners    = paste(mu_learners, collapse = ", "),
    g_learners     = paste(g_learners, collapse = ", "),
    kappa_learners = paste(kappa_learners, collapse = ", "),
    m_learners     = paste(m_learners, collapse = ", "),
    q_learners     = paste(q_learners, collapse = ", "),
    alpha = alpha,
    p1 = p1,
    p2 = p2,
    d1 = d1
  )

write.csv(resdf,
          file = "../output/sim_results.csv",
          row.names = FALSE)
