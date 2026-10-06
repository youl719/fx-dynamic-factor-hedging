# scripts/03_model_estimation.R

# 1. Load Required Packages
library(zoo)
library(xts)
library(dlm)

# 2. Load the Cleaned Data and LOO Data
env_data <- readRDS("data/01_cleaned_data.rds")
loo_env <- readRDS("data/02_loo_data.rds")

all_data <- env_data$all_data
train_window <- env_data$train_window
test_window <- env_data$test_window
currencies <- env_data$currencies

loo_data_list <- loo_env$loo_data_list
loo_static_betas_table <- loo_env$loo_static_betas_table

# 3. Define the Static OLS Baseline (Main Data)
static_betas <- sapply(currencies, function(ccy) {
  formula_str <- paste(ccy, "~ PC1_Dollar + PC2_Carry")
  model <- lm(as.formula(formula_str), data = as.data.frame(all_data[train_window]))
  return(coef(model))
})
static_betas_table <- t(static_betas)

# 4. Define the 252-day Rolling OLS Function
run_rolling_ols <- function(data_subset) {
  model <- lm(data_subset[, 1] ~ data_subset[, 2] + data_subset[, 3])
  return(coef(model))
}

get_rolling_oos <- function(ccy, data_source) {
  rolling_data <- data_source[, c(ccy, "PC1_Dollar", "PC2_Carry")]
  rolling_betas <- rollapply(
    data = rolling_data, width = 252, FUN = run_rolling_ols, 
    by.column = FALSE, align = "right"
  )
  colnames(rolling_betas) <- c("Alpha", "Beta_PC1", "Beta_PC2")
  # Return the lagged Out-of-Sample betas to prevent look-ahead bias
  return(lag.xts(as.xts(rolling_betas))[test_window])
}

# 5. Define the State-Space Kalman Filter Function
run_kalman_filter <- function(target_ccy, data_full, train_win, test_win) {
  Y_xts <- 100 * data_full[, target_ccy]
  X_xts <- data_full[, c("PC1_Dollar", "PC2_Carry")]
  
  train_Y <- as.matrix(coredata(Y_xts[train_win]))
  train_X <- as.matrix(coredata(X_xts[train_win, ]))
  full_Y <- as.matrix(coredata(Y_xts))
  full_X <- as.matrix(coredata(X_xts))
  
  # State-Space specification: Random Walk (Phi = I) with diffuse priors (C0)
  build_dlm <- function(par, predictors) {
    mod <- dlmModReg(predictors, addInt = TRUE,
                     dV = exp(par[1]),
                     dW = exp(par[2:4]))
    mod$C0 <- diag(1e7, 3) 
    return(mod)
  }
  
  build_dlm_train <- function(par) build_dlm(par, train_X)
  v0 <- var(as.numeric(train_Y), na.rm = TRUE)
  
  # Multi-start optimization to prevent local minima in Maximum Likelihood Estimation
  starts <- list(
    c(log(0.5 * v0), rep(log(1e-5), 3)),
    c(log(v0), rep(log(1e-4), 3)),
    c(log(0.1 * v0), rep(log(1e-6), 3))
  )
  
  fits <- lapply(starts, function(st) {
    try(dlmMLE(y = train_Y, parm = st, build = build_dlm_train, hessian = TRUE), silent = TRUE)
  })
  
  ok <- Filter(function(f) !inherits(f, "try-error") && f$convergence == 0, fits)
  
  if (length(ok) == 0) {
    best_fit <- dlmMLE(y = train_Y, parm = starts[[1]], build = build_dlm_train, hessian = TRUE)
  } else {
    best_fit <- ok[[which.min(sapply(ok, `[[`, "value"))]]
  }
  
  optimal_mod <- build_dlm(best_fit$par, full_X)
  filtered_out <- dlmFilter(full_Y, optimal_mod)
  
  states <- dropFirst(filtered_out$m) / 100
  states_xts <- as.xts(states, order.by = index(Y_xts))
  colnames(states_xts) <- c("Alpha", "Beta_PC1", "Beta_PC2")
  
  # Lag states for actionable out-of-sample hedging
  actionable_hedge_ratios <- lag.xts(states_xts)
  
  innov_xts <- xts(as.numeric(residuals(filtered_out, sd = FALSE)), order.by = index(Y_xts))
  innov_clean <- innov_xts[-(1:50)] # Drop initial diffuse prior volatility
  
  return(list(
    states = states_xts,
    oos_hedge = actionable_hedge_ratios[test_win],
    innov = innov_clean
  ))
}

# 6. Execute Models on Main Data
message("Running models on Main Factors...")
rolling_oos_all <- lapply(currencies, get_rolling_oos, data_source = all_data)
names(rolling_oos_all) <- currencies

kf_all <- lapply(currencies, function(ccy) {
  tryCatch(run_kalman_filter(ccy, all_data, train_window, test_window),
           error = function(e) { warning(ccy, ": ", conditionMessage(e)); NULL })
})
names(kf_all) <- currencies

# 7. Execute Models on LOO Data
message("Running models on LOO Factors...")
loo_rolling_oos_all <- lapply(currencies, function(ccy) get_rolling_oos(ccy, loo_data_list[[ccy]]))
names(loo_rolling_oos_all) <- currencies

loo_kf_all <- lapply(currencies, function(ccy) {
  tryCatch(run_kalman_filter(ccy, loo_data_list[[ccy]], train_window, test_window),
           error = function(e) { warning(ccy, ": ", conditionMessage(e)); NULL })
})
names(loo_kf_all) <- currencies

# 8. Save Model Fits
saveRDS(list(
  static_betas_table = static_betas_table,
  rolling_oos_all = rolling_oos_all,
  kf_all = kf_all,
  loo_rolling_oos_all = loo_rolling_oos_all,
  loo_kf_all = loo_kf_all
), file = "data/03_model_fits.rds")

print("Model estimation complete. Outputs saved to data/03_model_fits.rds")