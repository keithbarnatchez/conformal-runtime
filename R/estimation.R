# Propensity score trimming bounds
TRIM_LO <- 0.025
TRIM_HI <- 0.975

#' @title Weighted conformal prediction (internal)
#'
#' @description Sub-function called within est_r_dml to obtain initial r
#' estimate via weighted split conformal inference.
weighted_conformal <- function(Y,
                               A,
                               V, W,
                               S,
                               g_mod,
                               kappa_mod,
                               conf_score, 
                               score='abs_resid',
                               mu_learners = 'SL.glm.interaction',
                               alpha=0.1) {

 
  if (!is.null(W)) {
    X <- as.data.frame(cbind(V, W))
  } else {
    X <- as.data.frame(V)
  }
  
  n <- nrow(X)

  # Split conformal on the training/source data
  train_idx <- ifelse(S==0, 0, rbinom(n = length(S), size = 1, prob = 0.5))
  which_idx <- which(train_idx==1)

  # cal_idx should = 0 if S==0, and 1-train_idx if S==1
  cal_idx <- ifelse(S==0, 0, 1-train_idx)

  # Fit outcome in training data
  VA <- cbind(V,A)
  YXA <- cbind(Y,X,A)
  YVA <- cbind(Y,V,A)
  XA <- cbind(X,A)
  
  # Get g and kappa in calibration data (fit outside the function)
  g_hat <- SuperLearner::predict.SuperLearner(
    newdata = X,
    object = g_mod,
    onlySL=TRUE
  )$pred
  
  kappa_hat <- SuperLearner::predict.SuperLearner(
    newdata = V,
    object = kappa_mod,
    onlySL=TRUE
  )$pred
  
  # Trim propensity scores for stability
  g_hat <- pmin(pmax(g_hat, TRIM_LO), TRIM_HI)
  kappa_hat <- pmin(pmax(kappa_hat, TRIM_LO), TRIM_HI)

  # form normalized calibration weights
  of_weights <- (1/g_hat) * ((1-kappa_hat)/kappa_hat) / mean(cal_idx)
  of_weights_norm <-  n  * of_weights / sum(of_weights[cal_idx*A==1])
  max_score <- max(conf_score[cal_idx*A==1])
  
  # form the estimating eq to be solved
  # effectively finds the 1-alpha quantile of the weighted conf score dist
  regular_quantile <- function(r) {
    mean(conf_score <= r)
  }

  wgtd_quantile <- function(r) {
    mean(S * A * cal_idx * as.numeric(conf_score <= r) * of_weights_norm)
  }

  est_eq <- function(r) {
    return(
      wgtd_quantile(r) - (1-alpha)
    )
  }
  

  # Solve for r
  r <- min(uniroot(est_eq, c(min(conf_score), max(conf_score)))$root, max_score)
  
  return(r)

}

#' @title Estimate first stage CDF
#' 
#' @param r Conformity score cutoff
#' @param R Conformity scores
#' @param X Covariates for first stage regression
#' @param S Source indicator
#' @param A Treatment indicator
#' @param q_learners SuperLearner library for first stage regression
#' 
#' @return SuperLearner model for first stage regression
#' @export
est_q_a <- function(r,
                    R, X,
                    S,A,
                    q_learners = 'SL.glm.interaction') {

  ytilde <- as.numeric(R <= r)
  idx <- which(S==1)
  XA <- cbind(X,A)
  q_mod <- SuperLearner::SuperLearner(Y=ytilde[idx],
                                      X=XA[idx,,drop=FALSE],
                                      family = binomial(),
                                      SL.library = q_learners)

  return(q_mod)

}

#' @title Estimate second stage CDF
#'
#' @param q_a First stage conditional cdf
#' @param R Conformity scores
#' @param V Covariates for second stage regression
#' @param S Source indicator
#' @param A Treatment indicator
#' @param m_learners SuperLearner library for second stage regression
#'
#' @return SuperLearner model for second stage regression
#' @export
est_m_a <- function(q_a,
                    R, V,
                    S, A,
                    m_learners = 'SL.glm.interaction') {
  
  idx <- which(S==1)

  m_mod <- SuperLearner::SuperLearner(Y=q_a[idx],
                                      X=V[idx,,drop=FALSE],
                                      SL.library = m_learners)
  m_a = predict(m_mod)$pred
  
  return(m_mod)

}

#' @title Efficient influence function for the conformal CDF
#'
#' @param r Conformity score cutoff
#' @param R Conformity scores
#' @param A Treatment indicator
#' @param S Source indicator
#' @param W Runtime confounders (may be NULL)
#' @param V Always-observed covariates
#' @param g_a_hat Estimated propensity scores
#' @param kappa_hat Estimated source membership probabilities
#' @param m_a_hat Estimated second-stage CDF (E[q|V])
#' @param q_a_hat Estimated first-stage CDF
#' @param alpha Nominal miscoverage level
#'
#' @return Vector of influence function values
chi_a <- function(r,
                  R,A,S,W,V,
                  g_a_hat, kappa_hat,
                  m_a_hat,q_a_hat,
                  alpha) {

  chi_a_hat <- (1-S)*(m_a_hat - (1-alpha)) +
    S*(1-kappa_hat)/kappa_hat*(q_a_hat - m_a_hat) +
    (1-kappa_hat)*A*S/(kappa_hat*g_a_hat)*( as.numeric(R <= r) - q_a_hat )
  
  return(chi_a_hat)
}

#' @title Debiased conformal under runtime confounding
#' 
#'
#' @description Main function for performing debiased conformal prediction under
#' runtime confounding. 
#'
#' @param Y Vector containing the outcome variable
#' @param A Vector containing binary treatment values
#' @param V Dataframe containing covariates that are always available
#' @param W Dataframe containing covariates that are not available at runtime,
#' in the target data
#' @param S Vector with binary source population indicators: 0=target population,
#' 1=source population
#' @param score The conformal score function to use. Options are 'abs_resid' and 'quantile'
#' @param learners SuperLearner libraries to be used for prediction
#'
#' @return A dataframe with lower and upper bounds for Y(a) for each observation
#' in target population
#' @export
est_r_dml <- function(Y,
                      A,
                      V, W,
                      S,
                      score='abs_resid',
                      mu_learners = 'SL.glm.interaction',
                      g_learners = 'SL.glm.interaction',
                      kappa_learners = 'SL.glm.interaction',
                      m_learners = 'SL.glm.interaction',
                      q_learners = 'SL.glm.interaction',
                      alpha=0.1) {

  # Split data into disjoint sets D1 and D2
  n <- length(Y)
  D1 <- sample(1:n, size = ceiling(n/2), replace = FALSE)
  D2 <- setdiff(1:n, D1)
  Sidx <- which(S==1)
  D1S <- intersect(D1, Sidx)
  
  # Make training idx: all values of D1 where S also equals 1
  if (!is.null(W)) {
    X <- cbind(V,W)
  } else {
    X <- V
  }

  # Fit propensity scores on D1
  g_mod <- SuperLearner::SuperLearner(Y=A[D1S],
                                      X=X[D1S,,drop=FALSE],
                                      family = binomial(),
                                      SL.library = g_learners)

  g_hat <- SuperLearner::predict.SuperLearner(
    newdata = X,
    object = g_mod,
    onlySL=TRUE
  )$pred

  kappa_mod <- SuperLearner::SuperLearner(Y=S[D1],
                                          X=V[D1,,drop=FALSE],
                                          family=binomial(),
                                          SL.library = kappa_learners)
  kappa_hat <- SuperLearner::predict.SuperLearner(
    newdata = V,
    object = kappa_mod,
    onlySL=TRUE
  )$pred

  # Trim for stability
  g_hat <- pmin(pmax(g_hat, TRIM_LO), TRIM_HI)
  kappa_hat <- pmin(pmax(kappa_hat, TRIM_LO), TRIM_HI)
  
  VA <- cbind(V,A)
  XA <- cbind(X,A)
  YXA <- cbind(Y,X,A)
  YVA <- cbind(Y,V,A)
  
  if (score == 'abs_resid') {
    # Predict outcome on first fold and get scores
    mu_mod <- SuperLearner::SuperLearner(
      Y=Y[D1S],
      X=XA[D1S,,drop=FALSE],
      SL.library = mu_learners
    )
    
    # Extract predictions
    mu_hat <- SuperLearner::predict.SuperLearner(
      newdata = XA %>% mutate(A=1),
      object = mu_mod,
      onlySL=TRUE
    )$pred
    
    nu_hat <- mu_hat
    if (!is.null(W)) {
      # Further regress predictions on V
      nu_mod <- SuperLearner::SuperLearner(
        Y=mu_hat[D1S],
        X=V[D1S,,drop=FALSE],
        SL.library = mu_learners
      )
      
      nu_hat <- SuperLearner::predict.SuperLearner(
        newdata = V,
        object = nu_mod,
        onlySL=TRUE
      )$pred
    }
    
    R <- abs(Y-nu_hat)

  } 
  
  if (score == 'quantile') {
  
    a_weights <- A*S*(1-kappa_hat)/(kappa_hat*g_hat)
    
    if (!is.null(W)) {
      
    Q_mod <- ranger::ranger(
      Y ~ ., data=YVA[D1,],
      quantreg = TRUE,          
      case.weights = a_weights[D1],          
      num.trees = 500)
    } else {
      Q_mod <- ranger::ranger(
        Y ~ ., data=YVA[D1,],
        quantreg = TRUE,
        num.trees = 500)
    }
    
    Q_hat <- predict(Q_mod,
            data=YVA %>% mutate(A=1),
            type='quantiles',
            quantiles=c(alpha/2, 1-alpha/2))$predictions
    
    R <-  pmax(Q_hat[,1]-Y,Y-Q_hat[,2])
  }
  Rcal <- R[D1]
  max_R <- max(Rcal[A[D1]==1])
  
  # Further split first set into D11 and D12
  D11 <- sample(D1, size = length(D1)/2, replace = FALSE)
  D12 <- setdiff(D1, D11)

  # Perform weighted conformal on D11
  r_init <- weighted_conformal(Y[D11], A[D11], V[D11,,drop=FALSE],
                               W[D11,,drop=FALSE], S[D11],
                               g_mod, kappa_mod,
                               R[D11],
                               score, mu_learners,
                               alpha)

  
  # Train nuisance functions on D12
  # Estimate q and m
  q_mod <- est_q_a(r_init,
                   R[D12], X[D12,,drop=FALSE], S[D12], A[D12],
                   q_learners)
  q_hat <- SuperLearner::predict.SuperLearner(
    newdata = XA,
    object = q_mod,
    onlySL=TRUE
  )$pred

  if (!is.null(W)) {
  # Estimate m (E[q|V])
  m_mod <- est_m_a(q_hat[D12],
                   R = R[D12],
                   V = V[D12,,drop=FALSE],
                   S = S[D12],
                   A = A[D12],
                   m_learners)

  m_hat <- SuperLearner::predict.SuperLearner(
    newdata = V,
    object = m_mod,
    onlySL=TRUE
  )$pred
  } else {
    m_hat <- q_hat # will cancel out in EIF
  }

  # Form influence curve
  chi_hat <- chi_a(r_init,
                   R, A, S, W, V,
                   g_hat, kappa_hat,
                   m_hat, q_hat,
                   alpha)

  chi_hat_mean <- function(r) {
    mean(chi_a(r, R[D2], A[D2], S[D2], W[D2,], V[D2,],
               g_hat[D2], kappa_hat[D2], m_hat[D2], q_hat[D2], alpha))
  }

  # Find smallest r so that chi_hat_mean >= 0
  r_grid <- seq(min(R), max(R), length.out = 1e3)
  r_hat <- min(min(r_grid[which(sapply(r_grid, chi_hat_mean) >= 0)]), max_R)
  

  if (score == 'abs_resid') {
    # Now form prediction intervals for everyone
    int_L <- nu_hat - r_hat
    int_U <- nu_hat + r_hat
    
    int_L_init <- nu_hat - r_init
    int_U_init <- nu_hat + r_init
  }
  if (score == 'quantile') {
    # Now form prediction intervals for everyone
    int_L <- Q_hat[,1] - r_hat
    int_U <- Q_hat[,2] + r_hat
    
    int_L_init <- Q_hat[,1] - r_init
    int_U_init <- Q_hat[,2] + r_init
  }

  df <- cbind(X,
              data.frame(Y,A,S,
                         int_L, int_U,
                         int_L_init, int_U_init))
  
  if (score == 'abs_resid') {
    df <- df %>% mutate(nu_hat = nu_hat)
  }
  
  return(df)

}

#' @title Weighted conformal prediction with runtime confounding
#' 
#'
#' @description Weighted conformal prediction under runtime confounding.
#' Extends the method of Tibshirani et al. (2019) with importance weights
#' that account for distribution shift between source and target populations.
#'
#' @param Y Vector containing the outcome variable
#' @param A Vector containing binary treatment values
#' @param V Dataframe containing covariates that are always available
#' @param W Dataframe containing covariates that are not available at runtime,
#' in the target data
#' @param S Vector with binary source population indicators: 0=target population,
#' 1=source population
#' @param score The conformal score function to use. Options are 'abs_resid'
#' @param learners SuperLearner libraries to be used for prediction
#'
#' @return A dataframe with lower and upper bounds for E[Y(a)] for each observation
#' in target distribution
#' @export
wcp_rc <- function(Y,
                   A,
                   V, W,
                   S,
                   score='abs_resid',
                   mu_learners = 'SL.glm.interaction',
                   g_learners = 'SL.glm.interaction',
                   kappa_learners = 'SL.glm.interaction',
                   alpha=0.1) {
  
  n <- length(Y)
  if (!is.null(W)) {
    X <- cbind(V,W)
  } else {
    X <- V
  }

  # Get probability of treatment in source
  g_mod <- SuperLearner::SuperLearner(Y=A[S==1],
                                      X=X[S==1,,drop=FALSE],
                                      family = binomial(),
                                      SL.library = g_learners)
  
  # Get probability of source membership
  kappa_mod <- SuperLearner::SuperLearner(Y=S,
                                          X=V,
                                          family=binomial(),
                                          SL.library = kappa_learners)
  
  # Get g and kappa in calibration data
  g_hat <- SuperLearner::predict.SuperLearner(
    newdata = X,
    object = g_mod,
    onlySL=TRUE
  )$pred
  g_hat <- pmin(pmax(g_hat, TRIM_LO), TRIM_HI)

  kappa_hat <- SuperLearner::predict.SuperLearner(
    newdata = V,
    object = kappa_mod,
    onlySL=TRUE
  )$pred
  kappa_hat <- pmin(pmax(kappa_hat, TRIM_LO), TRIM_HI)
  
  # Split conformal on the training/source data
  # Training indices
  train_idx <- ifelse(S==0, 0, rbinom(n = length(S), size = 1, prob = 0.5))
  which_idx = which(train_idx==1)
  
  # Calibration indices
  cal_idx <- ifelse(S==0, 0, 1-train_idx)
  
  # Fit outcome in training data
  VA <- cbind(V,A)
  XA <- cbind(X,A)
  YVA <- cbind(Y,A,V)
  
  if (score == 'abs_resid') {
    
    # outcome regression on all X
    mu_mod <- SuperLearner::SuperLearner( 
      Y=Y[which_idx],
      X=XA[which_idx,, drop=FALSE],
      SL.library = mu_learners
    )
    
    # extract predictions
    mu_cal <- SuperLearner::predict.SuperLearner(
      newdata = XA %>% mutate(A=1),
      object = mu_mod,
      onlySL=TRUE
    )$pred
    
    nu_hat <- mu_cal 
    
    if (!is.null(W)) {
      # further regress on V
      nu_mod <- SuperLearner::SuperLearner(
        Y=mu_cal[which_idx],
        X=V[which_idx,, drop=FALSE],
        SL.library = mu_learners
      )
      
      # extract predictions
      nu_hat <- SuperLearner::predict.SuperLearner(
        newdata = V,
        object = nu_mod,
        onlySL=TRUE
      )$pred
    }
    conf_score <- abs(Y - nu_hat)
  }
  if (score == 'quantile') {

    a_weights <- A*S*(1-kappa_hat)/(kappa_hat*g_hat)
    
    if (is.null(W)) {
      Q_mod <- ranger::ranger(
        Y ~ ., data=YVA[which_idx,],
        quantreg = TRUE,
        num.trees = 500)
    } else {
      Q_mod <- ranger::ranger(
        Y ~ ., data=YVA[which_idx,],
        quantreg = TRUE,          
        case.weights = a_weights[which_idx],          
        num.trees = 500)
    }
    Q_hat <- predict(Q_mod,
                     data=YVA,
                     type='quantiles',
                     quantiles=c(alpha/2, 1-alpha/2))$predictions
    
    conf_score <- pmax(Q_hat[,1]-Y,Y-Q_hat[,2])
  }
  
  # form normalized calibration weights
  of_weights <- (1/g_hat) * ((1-kappa_hat)/kappa_hat) 
  of_weights_norm <-  n  * of_weights / sum(of_weights[cal_idx*A*S==1])
  
  # form the estimating eq to be solved
  # effectively finds the 1-alpha quantile of the weighted conf score dist
  regular_quantile <- function(r) {
    mean(conf_score <= r)
  }
  
  wgtd_quantile <- function(r) {
    mean(S * A * cal_idx * as.numeric(conf_score <= r) * of_weights_norm)
  }
  
  est_eq <- function(r) {
    return(
      wgtd_quantile(r) - (1-alpha)
    )
  }
  
  # Collect scores and weights into a df
  score_df <- data.frame(scores=conf_score,
                         g_hat=g_hat,
                         kappa_hat=kappa_hat,
                         weight = (1/g_hat)*(1-kappa_hat)/kappa_hat,
                         cal_idx=cal_idx)
  score_df <- score_df %>% filter(cal_idx==1,A==1)
  
  get_rval <- function(new_obs,
                       alpha) {
    
    data_scores <- score_df$scores
    norm_weights <- score_df$weight / c(sum(score_df$weight) + new_obs$weight)
    new_weight <- new_obs$weight / c(sum(score_df$weight) + new_obs$weight)
    scores <- c(data_scores,max(score_df$scores)) 
    norm_weights <- c(norm_weights,new_weight)
    
    rval <- wgtd_quantile_f(scores,1-alpha,norm_weights)

    return(rval)
    
  }

  # Solve for r
  the_idx = which(A*S*cal_idx==1)
  
  df <- cbind(X,
              data.frame(Y,A,S))
  df$weight <- (1-kappa_hat)/kappa_hat
  df$kappa_hat <- kappa_hat
  df <- df %>% filter(S==0) 
  df$rvals <- rep(NA,nrow(df))

  
  for (rr in 1:nrow(df)) {
    new_obs <- df[rr,]
    rval <- get_rval(new_obs,alpha)
    df[rr,'rvals'] <- rval
  }
  
  if (score == 'abs_resid') {
    df$int_L <- nu_hat[(1-S)==1] - df$rvals
    df$int_U <- nu_hat[(1-S)==1] + df$rvals
  }
  if (score == 'quantile') {
    df$int_L <- Q_hat[(1-S)==1,1] - df$rvals
    df$int_U <- Q_hat[(1-S)==1,2] + df$rvals
  }
  
  return(df)
  
}

#' @title Weighted quantile function
#'
#' @param v Numeric vector of values
#' @param prob Desired quantile probability
#' @param weights Numeric vector of weights (same length as v)
#'
#' @return The weighted quantile value
wgtd_quantile_f <- function(v, prob, weights) {
  
  # Sort the values and weights
  sorted_indices <- order(v)
  sorted_values <- v[sorted_indices]
  sorted_weights <- weights[sorted_indices]
  
  # Calculate the cumulative sum of weights
  cumulative_weights <- cumsum(sorted_weights)
  
  # Find the smallest index where the cumulative weight exceeds the desired probability
  target_index <- which(cumulative_weights >= prob * sum(sorted_weights))[1]
  
  # Return the corresponding value
  return(sorted_values[target_index])
  
}

