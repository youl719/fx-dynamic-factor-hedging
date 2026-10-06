# 1. Load Required Packages
library(quantmod)
library(zoo)
library(dlm)
library(PerformanceAnalytics)

# 2. Define FRED Series Tickers for G10 Currencies
tickers <- c(
  EUR = "DEXUSEU",
  GBP = "DEXUSUK",
  AUD = "DEXUSAL",
  NZD = "DEXUSNZ",
  JPY = "DEXJPUS", # Inverted (JPY per USD)
  CAD = "DEXCAUS", # Inverted (CAD per USD)
  CHF = "DEXSZUS", # Inverted (CHF per USD)
  SEK = "DEXSDUS", # Inverted (SEK per USD)
  NOK = "DEXNOUS"  # Inverted (NOK per USD)
)

# 3. Download Data via quantmod
getSymbols(tickers, src = "FRED", env = globalenv())

# Merge into a single xts (time series) object
fx_prices <- merge(DEXUSEU, DEXUSUK, DEXUSAL, DEXUSNZ, DEXJPUS, DEXCAUS, DEXSZUS, DEXSDUS, DEXNOUS)
colnames(fx_prices) <- names(tickers)

# 4. Clean Missing Data (Fixed Holiday Handling)
fx_prices <- fx_prices[rowSums(!is.na(fx_prices)) > 0, ]
fx_prices <- na.locf(fx_prices)

# 5. Standardize Quoting & Calculate Log Returns
fx_prices$JPY <- 1 / fx_prices$JPY
fx_prices$CAD <- 1 / fx_prices$CAD
fx_prices$CHF <- 1 / fx_prices$CHF
fx_prices$SEK <- 1 / fx_prices$SEK
fx_prices$NOK <- 1 / fx_prices$NOK

fx_returns <- diff(log(fx_prices))
fx_returns <- na.omit(fx_returns)

# 6. Split into In-Sample (Train) and Out-of-Sample (Test)
train_window <- "2005-01-01/2018-12-31"
test_window <- "2019-01-01/2025-12-31"

train_returns <- fx_returns[train_window]
test_returns <- fx_returns[test_window]

# 7. Pure PCA on the Training Data (Fixed Sign Indeterminacy)
pca_train <- prcomp(train_returns, center = TRUE, scale. = TRUE)

rot <- pca_train$rotation[, 1:2]
if (mean(rot[, 1]) < 0) rot[, 1] <- -rot[, 1]
if (mean(rot[c("AUD", "NZD"), 2]) < 0) rot[, 2] <- -rot[, 2]

# Verification: Print Variance Explained and Sign-Corrected Loadings
print("PCA Variance Explained:")
print(summary(pca_train)$importance[, 1:2])

print("PCA Loadings (Rotation Matrix):")
print(round(rot, 3))

scale_fx <- function(x) scale(x, center = pca_train$center, scale = pca_train$scale)

factors_train <- as.xts(scale_fx(train_returns) %*% rot, order.by = index(train_returns))
colnames(factors_train) <- c("PC1_Dollar", "PC2_Carry")

# 8. Out-of-Sample Projection
test_returns_scaled <- scale_fx(test_returns)
stopifnot(!anyNA(test_returns_scaled))

factors_test <- as.xts(test_returns_scaled %*% rot, order.by = index(test_returns))
colnames(factors_test) <- c("PC1_Dollar", "PC2_Carry")

# Verification: Plot the COVID-19 Liquidity Shock
plot.zoo(factors_test["2020-01-01/2020-06-30"]$PC1_Dollar,
         main = "OOS PC1 (Dollar Factor) during COVID-19",
         ylab = "PC1",
         col = "blue")

# 9. Combine Train and Test for Continuous Modeling
all_factors <- rbind(factors_train, factors_test)
all_returns <- rbind(train_returns, test_returns)
all_data <- na.omit(merge(all_returns, all_factors))
train_data <- all_data[train_window]

# 10. Static OLS: Full Sample (2005-2018)
currencies <- c("EUR", "GBP", "AUD", "NZD", "JPY", "CAD", "CHF", "SEK", "NOK")

static_betas <- sapply(currencies, function(ccy) {
  formula_str <- paste(ccy, "~ PC1_Dollar + PC2_Carry")
  model <- lm(as.formula(formula_str), data = as.data.frame(train_data))
  return(coef(model))
})

static_betas_table <- t(static_betas)

# 11. Rolling-Window OLS: 252 Days (Example: AUD)
run_rolling_ols <- function(data_subset) {
  model <- lm(data_subset[, 1] ~ data_subset[, 2] + data_subset[, 3])
  return(coef(model))
}

aud_rolling_data <- all_data[, c("AUD", "PC1_Dollar", "PC2_Carry")]

rolling_betas_aud <- rollapply(
  data = aud_rolling_data,
  width = 252,
  FUN = run_rolling_ols,
  by.column = FALSE,
  align = "right"
)

colnames(rolling_betas_aud) <- c("Alpha", "Beta_PC1", "Beta_PC2")
rolling_betas_aud <- as.xts(rolling_betas_aud)

# Extract lagged OOS period
oos_rolling_betas <- lag.xts(rolling_betas_aud)[test_window]

# Verification: Hard Assertion for Look-Ahead Bias
first_test_day <- index(oos_rolling_betas)[1]
last_train_day <- tail(index(rolling_betas_aud[train_window]), 1)
stopifnot(all.equal(as.numeric(oos_rolling_betas[first_test_day]),
                    as.numeric(rolling_betas_aud[last_train_day])))

# 12. Visualizing the Instability
plot.zoo(rolling_betas_aud$Beta_PC1[train_window],
         type = "l",
         col = "darkred",
         main = "AUD In-Sample Rolling Beta vs. Static Beta",
         ylab = "PC1 Beta",
         xlab = "Date")
abline(h = static_betas_table["AUD", "PC1_Dollar"], col = "blue", lwd = 2, lty = 2)

# 13. State-Space Kalman Filter Implementation (Final Validation Polish)
run_kalman_filter <- function(target_ccy, data_full, train_win, test_win) {
  Y_xts <- 100 * data_full[, target_ccy]
  X_xts <- data_full[, c("PC1_Dollar", "PC2_Carry")]
  
  train_Y <- as.matrix(coredata(Y_xts[train_win]))
  train_X <- as.matrix(coredata(X_xts[train_win, ]))
  full_Y <- as.matrix(coredata(Y_xts))
  full_X <- as.matrix(coredata(X_xts))
  
  build_dlm <- function(par, predictors) {
    mod <- dlmModReg(predictors, addInt = TRUE,
                     dV = exp(par[1]),
                     dW = exp(par[2:4]))
    mod$C0 <- diag(1e7, 3)
    return(mod)
  }
  
  build_dlm_train <- function(par) build_dlm(par, train_X)
  
  v0 <- var(as.numeric(train_Y), na.rm = TRUE)
  
  starts <- list(
    c(log(0.5 * v0), rep(log(1e-5), 3)),
    c(log(v0), rep(log(1e-4), 3)),
    c(log(0.1 * v0), rep(log(1e-6), 3))
  )
  
  # Evaluate all starting parameters
  fits <- lapply(starts, function(st) {
    try(dlmMLE(y = train_Y, parm = st, build = build_dlm_train, hessian = TRUE), silent = TRUE)
  })
  
  ok <- Filter(function(f) !inherits(f, "try-error") && f$convergence == 0, fits)
  
  if (length(ok) == 0) {
    warning(paste("All MLE optimizer starts failed or did not converge for", target_ccy))
    best_fit <- dlmMLE(y = train_Y, parm = starts[[1]], build = build_dlm_train, hessian = TRUE)
    all_fits_summary <- NULL
  } else {
    all_fits_summary <- sapply(ok, function(f) c(value = f$value, exp(f$par)))
    best_fit <- ok[[which.min(sapply(ok, `[[`, "value"))]]
  }
  
  optimal_mod <- build_dlm(best_fit$par, full_X)
  filtered_out <- dlmFilter(full_Y, optimal_mod)
  
  states <- dropFirst(filtered_out$m) / 100
  states_xts <- as.xts(states, order.by = index(Y_xts))
  colnames(states_xts) <- c("Alpha", "Beta_PC1", "Beta_PC2")
  
  actionable_hedge_ratios <- lag.xts(states_xts)
  
  innov_xts <- xts(as.numeric(residuals(filtered_out, sd = FALSE)), order.by = index(Y_xts))
  innov_clean <- innov_xts[-(1:50)]
  
  return(list(
    mle_fit = best_fit,
    all_starts = all_fits_summary,
    states = states_xts,
    oos_hedge = actionable_hedge_ratios[test_win],
    innov = innov_clean
  ))
}

# 14. Execute Kalman Filter
kf_aud <- run_kalman_filter("AUD", all_data, train_window, test_window)

# 15. Sanity Checks & Diagnostics
print("--- Multi-Start Optimization Comparison [value, V, W_alpha, W_PC1, W_PC2] ---")
print(kf_aud$all_starts)

print("--- Standard Errors on Log Scale ---")
tryCatch({
  inv_hessian <- solve(kf_aud$mle_fit$hessian)
  diag_inv <- diag(inv_hessian)
  if (any(diag_inv < 0)) {
    print("Hessian not positive definite - negative values on diagonal (NaN SEs)")
  } else {
    print(sqrt(diag_inv))
  }
}, error = function(e) print("Hessian singular - variance likely at boundary"))

print("--- End-of-Training Kalman Betas ---")
print(tail(kf_aud$states[train_window], 1))

print("--- Static OLS Betas (For Magnitude Comparison) ---")
print(static_betas_table["AUD", ])

print("--- Summary: Kalman vs Rolling OOS Hedge Ratios (COVID-19) ---")
k <- kf_aud$oos_hedge$Beta_PC1["2020-01-01/2020-06-30"]
r <- oos_rolling_betas$Beta_PC1["2020-01-01/2020-06-30"]
print(summary(cbind(Kalman = k, Rolling = r)))

# 16. Verify Out-of-Sample COVID-19 Liquidity Shock Adaptability
plot.zoo(as.zoo(k), ylim = range(c(k, r), na.rm = TRUE), col = "blue", lwd = 2,
         main = "AUD PC1 Beta: Kalman vs 252d Rolling OLS (COVID-19)",
         ylab = "PC1 Hedge Ratio", xlab = "Date")
lines(as.zoo(r), col = "darkred", lwd = 2)
legend("bottomright", legend = c("Kalman Filter (m_t)", "252-day Rolling OLS"),
       col = c("blue", "darkred"), lwd = 2,bg="transparent")

#---
# 17. Multi-Asset Execution: Rolling OLS & Kalman Filter
library(PerformanceAnalytics)

currencies <- c("EUR", "GBP", "AUD", "NZD", "JPY", "CAD", "CHF", "SEK", "NOK")

get_rolling_oos <- function(ccy) {
  rolling_data <- all_data[, c(ccy, "PC1_Dollar", "PC2_Carry")]
  rolling_betas <- rollapply(
    data = rolling_data, width = 252, FUN = run_rolling_ols, 
    by.column = FALSE, align = "right"
  )
  colnames(rolling_betas) <- c("Alpha", "Beta_PC1", "Beta_PC2")
  return(lag.xts(as.xts(rolling_betas))[test_window])
}
rolling_oos_all <- lapply(currencies, get_rolling_oos)
names(rolling_oos_all) <- currencies

kf_all <- lapply(currencies, function(ccy) {
  message("Kalman: ", ccy)
  tryCatch(
    run_kalman_filter(ccy, all_data, train_window, test_window),
    error = function(e) { warning(ccy, ": ", conditionMessage(e)); NULL }
  )
})
names(kf_all) <- currencies


# 18. Model Diagnostics: Ljung-Box Test (Dynamic Whitening & Heteroskedasticity)
lb_pval <- function(x, lag = 10) Box.test(as.numeric(x), lag = lag, type = "Ljung-Box")$p.value

lb_results <- sapply(currencies, function(ccy) {
  if (is.null(kf_all[[ccy]])) {
    return(setNames(rep(NA_real_, 6), 
                    c("Static_IS", "Kalman_IS", "Static_OOS", "Kalman_OOS", "Static_OOS_Sq", "Kalman_OOS_Sq")))
  }
  
  mod <- lm(as.formula(paste(ccy, "~ PC1_Dollar + PC2_Carry")), data = as.data.frame(train_data))
  b <- coef(mod)
  d <- all_data[test_window]
  static_oos <- as.numeric(d[, ccy]) - b[1] - b[2] * as.numeric(d$PC1_Dollar) - b[3] * as.numeric(d$PC2_Carry)
  
  innov <- kf_all[[ccy]]$innov
  
  c(Static_IS     = lb_pval(residuals(mod)),
    Kalman_IS     = lb_pval(innov[train_window]),
    Static_OOS    = lb_pval(static_oos),
    Kalman_OOS    = lb_pval(innov[test_window]),
    Static_OOS_Sq = lb_pval(static_oos^2),
    Kalman_OOS_Sq = lb_pval(innov[test_window]^2))
})
print("--- Ljung-Box Test P-Values (Lag 10) ---")
print(round(t(lb_results), 4))


# 19. Out-of-Sample Hedged Returns Construction (Alpha Excluded)
oos_factors <- all_data[test_window, c("PC1_Dollar", "PC2_Carry")]

calculate_hedged_returns <- function(ccy) {
  if (is.null(kf_all[[ccy]])) return(NULL)
  
  y_oos <- all_data[test_window, ccy]
  
  # Strict factor hedge: excludes local alpha drift
  hedge <- function(y, b1, b2) {
    y - (b1 * oos_factors$PC1_Dollar + b2 * oos_factors$PC2_Carry)
  }
  
  unhedged <- y_oos
  
  static_b <- static_betas_table[ccy, ]
  static_hedged <- hedge(y_oos, static_b["PC1_Dollar"], static_b["PC2_Carry"])
  
  roll_b <- rolling_oos_all[[ccy]]
  roll_hedged <- hedge(y_oos, roll_b$Beta_PC1, roll_b$Beta_PC2)
  
  kf_b <- kf_all[[ccy]]$oos_hedge
  kf_hedged <- hedge(y_oos, kf_b$Beta_PC1, kf_b$Beta_PC2)
  
  out <- merge(unhedged, static_hedged, roll_hedged, kf_hedged)
  colnames(out) <- c("Unhedged", "Static", "Rolling", "Kalman")
  
  return(na.omit(out))
}

oos_hedged_returns <- lapply(currencies, calculate_hedged_returns)
names(oos_hedged_returns) <- currencies


# 20. Final Evaluation Metrics & Master Tables
evaluate_performance <- function(ccy) {
  ret_data <- oos_hedged_returns[[ccy]]
  if (is.null(ret_data)) return(NULL)
  
  # Full OOS Window: Variance Reduction & Empirical VaR
  v <- apply(ret_data, 2, var)
  var_reduction_full <- 100 * (1 - v / v["Unhedged"])
  
  var_95_full <- -100 * apply(ret_data, 2, quantile, probs = 0.05)
  var_99_full <- -100 * apply(ret_data, 2, quantile, probs = 0.01)
  
  # Feb-May 2020 Specific Stress Period
  covid_window <- "2020-02-01/2020-05-31"
  covid_ret <- ret_data[covid_window]
  
  v_c <- apply(covid_ret, 2, var)
  var_red_covid <- 100 * (1 - v_c / v_c["Unhedged"])
  
  sse_covid <- colSums(covid_ret^2)
  sse_ratio_covid <- 100 * sse_covid / sse_covid["Unhedged"]
  
  # Extract names to preserve xts objects during iteration
  max_dd_covid <- 100 * sapply(colnames(covid_ret), function(nm) {
    maxDrawdown(covid_ret[, nm], geometric = FALSE)
  })
  
  metrics <- rbind(Variance_Reduction_Full = var_reduction_full,
                   VaR_95_Full = var_95_full,
                   VaR_99_Full = var_99_full,
                   Variance_Reduction_COVID = var_red_covid,
                   SSE_Ratio_COVID = sse_ratio_covid,
                   Max_Drawdown_COVID = max_dd_covid)
  
  metrics["Variance_Reduction_Full", "Unhedged"] <- 0 
  metrics["Variance_Reduction_COVID", "Unhedged"] <- 0 
  
  return(round(metrics, 2))
}

eval_results <- lapply(currencies, evaluate_performance)
names(eval_results) <- currencies

# Compile Master Tables
valid_ccys <- currencies[!sapply(eval_results, is.null)]

cat("\n--- FULL OOS: Variance Reduction (%) ---\n")
print(t(sapply(valid_ccys, function(ccy) {
  eval_results[[ccy]]["Variance_Reduction_Full", c("Static", "Rolling", "Kalman")]
})))

cat("\n--- COVID-19 SHOCK (Feb-May 2020): Variance Reduction (%) ---\n")
print(t(sapply(valid_ccys, function(ccy) {
  eval_results[[ccy]]["Variance_Reduction_COVID", c("Static", "Rolling", "Kalman")]
})))

cat("\n--- COVID-19 SHOCK (Feb-May 2020): SSE Ratio (% of Unhedged Remaining) ---\n")
print(t(sapply(valid_ccys, function(ccy) {
  eval_results[[ccy]]["SSE_Ratio_COVID", c("Static", "Rolling", "Kalman")]
})))

cat("\n--- COVID-19 SHOCK (Feb-May 2020): Maximum Drawdown (%) ---\n")
print(t(sapply(valid_ccys, function(ccy) {
  eval_results[[ccy]]["Max_Drawdown_COVID", c("Unhedged", "Static", "Rolling", "Kalman")]
})))

cat("\n--- FULL OOS (2019-2025): 95% Historical VaR (%) ---\n")
print(t(sapply(valid_ccys, function(ccy) {
  eval_results[[ccy]]["VaR_95_Full", c("Unhedged", "Static", "Rolling", "Kalman")]
})))

cat("\n--- FULL OOS (2019-2025): 99% Historical VaR (%) ---\n")
print(t(sapply(valid_ccys, function(ccy) {
  eval_results[[ccy]]["VaR_99_Full", c("Unhedged", "Static", "Rolling", "Kalman")]
})))

#---
# 21. Leave-One-Out (LOO) PCA, Data Structuring, & Correlation Alignment
build_loo <- function(ccy) {
  others <- setdiff(currencies, ccy)
  pca <- prcomp(train_returns[, others], center = TRUE, scale. = TRUE)
  rot <- pca$rotation[, 1:2]
  
  sc <- function(x) scale(x[, others], center = pca$center, scale = pca$scale)
  
  f_tr <- as.xts(sc(train_returns) %*% rot, order.by = index(train_returns))
  f_te <- as.xts(sc(test_returns) %*% rot, order.by = index(test_returns))
  
  f <- rbind(f_tr, f_te)
  colnames(f) <- c("PC1_Dollar", "PC2_Carry")
  
  # Align signs to the main-model factors (same training dates)
  main_tr <- all_data[train_window, c("PC1_Dollar", "PC2_Carry")]
  stopifnot(nrow(main_tr) == nrow(f[train_window]))
  
  cors <- c(PC1 = cor(as.numeric(f[train_window][, 1]), as.numeric(main_tr[, 1])),
            PC2 = cor(as.numeric(f[train_window][, 2]), as.numeric(main_tr[, 2])))
  
  f[, 1] <- sign(cors[1]) * f[, 1]
  f[, 2] <- sign(cors[2]) * f[, 2]
  
  dat <- na.omit(merge(rbind(train_returns[, ccy], test_returns[, ccy]), f))
  colnames(dat) <- c(ccy, "PC1_Dollar", "PC2_Carry")
  stopifnot(nrow(dat) == nrow(all_data))
  
  mod <- lm(as.formula(paste(ccy, "~ PC1_Dollar + PC2_Carry")),
            data = as.data.frame(dat[train_window]))
  
  list(data = dat, static = coef(mod), cors = abs(cors))
}

loo_objs <- setNames(lapply(currencies, build_loo), currencies)
loo_data_list <- lapply(loo_objs, `[[`, "data")
loo_static_betas_table <- t(sapply(loo_objs, `[[`, "static"))

cat("\n--- LOO Factor Similarity (Absolute Correlation with Main Factors) ---\n")
print(round(t(sapply(loo_objs, `[[`, "cors")), 3))


# 22. LOO Multi-Asset Execution: Rolling OLS & Kalman
loo_rolling_oos_all <- lapply(currencies, function(ccy) {
  rolling_data <- loo_data_list[[ccy]][, c(ccy, "PC1_Dollar", "PC2_Carry")]
  rolling_betas <- rollapply(
    data = rolling_data, width = 252, FUN = run_rolling_ols, 
    by.column = FALSE, align = "right"
  )
  colnames(rolling_betas) <- c("Alpha", "Beta_PC1", "Beta_PC2")
  return(lag.xts(as.xts(rolling_betas))[test_window])
})
names(loo_rolling_oos_all) <- currencies

loo_kf_all <- lapply(currencies, function(ccy) {
  message("LOO Kalman: ", ccy)
  tryCatch(
    run_kalman_filter(ccy, loo_data_list[[ccy]], train_window, test_window),
    error = function(e) { warning(ccy, ": ", conditionMessage(e)); NULL }
  )
})
names(loo_kf_all) <- currencies


# 23. LOO Out-of-Sample Hedged Returns Construction
calc_loo_hedged <- function(ccy) {
  if (is.null(loo_kf_all[[ccy]])) return(NULL)
  
  d_test <- loo_data_list[[ccy]][test_window]
  y_oos <- d_test[, ccy]
  oos_factors <- d_test[, c("PC1_Dollar", "PC2_Carry")]
  
  hedge <- function(y, b1, b2) {
    y - (b1 * oos_factors$PC1_Dollar + b2 * oos_factors$PC2_Carry)
  }
  
  unhedged <- y_oos
  static_b <- loo_static_betas_table[ccy, ]
  
  static_hedged <- hedge(y_oos, static_b["PC1_Dollar"], static_b["PC2_Carry"])
  roll_hedged <- hedge(y_oos, loo_rolling_oos_all[[ccy]]$Beta_PC1, loo_rolling_oos_all[[ccy]]$Beta_PC2)
  kf_hedged <- hedge(y_oos, loo_kf_all[[ccy]]$oos_hedge$Beta_PC1, loo_kf_all[[ccy]]$oos_hedge$Beta_PC2)
  
  out <- merge(unhedged, static_hedged, roll_hedged, kf_hedged)
  colnames(out) <- c("Unhedged", "Static", "Rolling", "Kalman")
  
  return(na.omit(out))
}

loo_hedged_returns <- lapply(currencies, calc_loo_hedged)
names(loo_hedged_returns) <- currencies


# 24. LOO Final Evaluation Metrics & Scale-Free Gap Analysis
loo_eval_results <- lapply(currencies, function(ccy) {
  ret_data <- loo_hedged_returns[[ccy]]
  if (is.null(ret_data)) return(NULL)
  
  v <- apply(ret_data, 2, var)
  var_reduction_full <- 100 * (1 - v / v["Unhedged"])
  
  covid_window <- "2020-02-01/2020-05-31"
  covid_ret <- ret_data[covid_window]
  v_c <- apply(covid_ret, 2, var)
  var_red_covid <- 100 * (1 - v_c / v_c["Unhedged"])
  
  metrics <- rbind(Variance_Reduction_Full = var_reduction_full,
                   Variance_Reduction_COVID = var_red_covid)
  
  metrics["Variance_Reduction_Full", "Unhedged"] <- 0 
  metrics["Variance_Reduction_COVID", "Unhedged"] <- 0 
  
  return(round(metrics, 2))
})
names(loo_eval_results) <- currencies

loo_valid_ccys <- currencies[!sapply(loo_eval_results, is.null)]

# Calculate Scale-Free Gap: log(Var_Rolling / Var_Kalman)
# Positive values indicate Var_Rolling > Var_Kalman (Kalman outperforms)
get_log_var_ratio <- function(ret_list, valid_c, window = NULL) {
  sapply(valid_c, function(c) {
    dat <- ret_list[[c]]
    if (!is.null(window)) dat <- dat[window]
    v <- apply(dat, 2, var)
    log(v["Rolling"] / v["Kalman"])
  })
}

log_ratio_full <- rbind(
  Main = get_log_var_ratio(oos_hedged_returns, valid_ccys),
  LOO  = get_log_var_ratio(loo_hedged_returns, loo_valid_ccys)
)

log_ratio_covid <- rbind(
  Main = get_log_var_ratio(oos_hedged_returns, valid_ccys, window = "2020-02-01/2020-05-31"),
  LOO  = get_log_var_ratio(loo_hedged_returns, loo_valid_ccys, window = "2020-02-01/2020-05-31")
)

cat("\n=== LEAVE-ONE-OUT (LOO) ROBUSTNESS CHECK ===\n")

cat("\n--- LOO FULL OOS: Variance Reduction (%) ---\n")
print(t(sapply(loo_valid_ccys, function(ccy) loo_eval_results[[ccy]]["Variance_Reduction_Full", c("Static", "Rolling", "Kalman")])))

cat("\n--- LOO COVID-19 SHOCK (Feb-May 2020): Variance Reduction (%) ---\n")
print(t(sapply(loo_valid_ccys, function(ccy) loo_eval_results[[ccy]]["Variance_Reduction_COVID", c("Static", "Rolling", "Kalman")])))

cat("\n--- SCALE-FREE GAP ANALYSIS: log(Var_Rolling / Var_Kalman) ---\n")
cat("    (Values > 0 indicate Kalman outperforms Rolling)\n")
cat("\n--- Full OOS Window ---\n")
print(round(log_ratio_full, 4))

cat("\n--- COVID-19 Shock Window ---\n")
print(round(log_ratio_covid, 4))