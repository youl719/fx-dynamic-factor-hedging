# scripts/04_evaluation_metrics.R

# 1. Load Required Packages
library(zoo)
library(xts)
library(PerformanceAnalytics)

# 2. Load the Data and Model Fits
env_data <- readRDS("data/01_cleaned_data.rds")
loo_env <- readRDS("data/02_loo_data.rds")
fits <- readRDS("data/03_model_fits.rds")

all_data <- env_data$all_data
train_data <- all_data[env_data$train_window]
train_window <- env_data$train_window
test_window <- env_data$test_window
currencies <- env_data$currencies

loo_data_list <- loo_env$loo_data_list
static_betas_table <- fits$static_betas_table
loo_static_betas_table <- loo_env$loo_static_betas_table

# 3. Model Diagnostics: Ljung-Box Test (Heteroskedasticity & Whitening)
lb_pval <- function(x, lag = 10) Box.test(as.numeric(x), lag = lag, type = "Ljung-Box")$p.value

lb_results <- sapply(currencies, function(ccy) {
  if (is.null(fits$kf_all[[ccy]])) return(rep(NA_real_, 6))
  
  mod <- lm(as.formula(paste(ccy, "~ PC1_Dollar + PC2_Carry")), data = as.data.frame(train_data))
  b <- coef(mod)
  d <- all_data[test_window]
  static_oos <- as.numeric(d[, ccy]) - b[1] - b[2] * as.numeric(d$PC1_Dollar) - b[3] * as.numeric(d$PC2_Carry)
  
  innov <- fits$kf_all[[ccy]]$innov
  
  c(Static_IS     = lb_pval(residuals(mod)),
    Kalman_IS     = lb_pval(innov[train_window]),
    Static_OOS    = lb_pval(static_oos),
    Kalman_OOS    = lb_pval(innov[test_window]),
    Static_OOS_Sq = lb_pval(static_oos^2),
    Kalman_OOS_Sq = lb_pval(innov[test_window]^2))
})
rownames(lb_results) <- c("Static_IS", "Kalman_IS", "Static_OOS", "Kalman_OOS", "Static_OOS_Sq", "Kalman_OOS_Sq")

# 4. Out-of-Sample Hedged Returns Construction (Main & LOO)
oos_factors <- all_data[test_window, c("PC1_Dollar", "PC2_Carry")]

calc_hedged <- function(ccy, is_loo = FALSE) {
  if (!is_loo && is.null(fits$kf_all[[ccy]])) return(NULL)
  if (is_loo && is.null(fits$loo_kf_all[[ccy]])) return(NULL)
  
  if (is_loo) {
    y_oos <- loo_data_list[[ccy]][test_window, ccy]
    fac <- loo_data_list[[ccy]][test_window, c("PC1_Dollar", "PC2_Carry")]
    stat_b <- loo_static_betas_table[ccy, ]
    roll_b <- fits$loo_rolling_oos_all[[ccy]]
    kf_b <- fits$loo_kf_all[[ccy]]$oos_hedge
  } else {
    y_oos <- all_data[test_window, ccy]
    fac <- oos_factors
    stat_b <- static_betas_table[ccy, ]
    roll_b <- fits$rolling_oos_all[[ccy]]
    kf_b <- fits$kf_all[[ccy]]$oos_hedge
  }
  
  hedge <- function(y, b1, b2) y - (b1 * fac$PC1_Dollar + b2 * fac$PC2_Carry)
  
  out <- merge(y_oos, 
               hedge(y_oos, stat_b["PC1_Dollar"], stat_b["PC2_Carry"]),
               hedge(y_oos, roll_b$Beta_PC1, roll_b$Beta_PC2),
               hedge(y_oos, kf_b$Beta_PC1, kf_b$Beta_PC2))
  colnames(out) <- c("Unhedged", "Static", "Rolling", "Kalman")
  return(na.omit(out))
}

hedged_returns_main <- lapply(currencies, calc_hedged, is_loo = FALSE)
names(hedged_returns_main) <- currencies

hedged_returns_loo <- lapply(currencies, calc_hedged, is_loo = TRUE)
names(hedged_returns_loo) <- currencies

# 5. Calculate Final Evaluation Metrics
evaluate_performance <- function(ret_data) {
  if (is.null(ret_data)) return(NULL)
  
  v <- apply(ret_data, 2, var)
  var_reduction_full <- 100 * (1 - v / v["Unhedged"])
  var_95_full <- -100 * apply(ret_data, 2, quantile, probs = 0.05)
  var_99_full <- -100 * apply(ret_data, 2, quantile, probs = 0.01)
  
  covid_window <- "2020-02-01/2020-05-31"
  covid_ret <- ret_data[covid_window]
  
  v_c <- apply(covid_ret, 2, var)
  var_red_covid <- 100 * (1 - v_c / v_c["Unhedged"])
  
  sse_covid <- colSums(covid_ret^2)
  sse_ratio_covid <- 100 * sse_covid / sse_covid["Unhedged"]
  
  max_dd_covid <- 100 * sapply(colnames(covid_ret), function(nm) maxDrawdown(covid_ret[, nm], geometric = FALSE))
  
  metrics <- rbind(Variance_Reduction_Full = var_reduction_full,
                   VaR_95_Full = var_95_full,
                   VaR_99_Full = var_99_full,
                   Variance_Reduction_COVID = var_red_covid,
                   SSE_Ratio_COVID = sse_ratio_covid,
                   Max_Drawdown_COVID = max_dd_covid)
  metrics[c("Variance_Reduction_Full", "Variance_Reduction_COVID"), "Unhedged"] <- 0 
  return(round(metrics, 2))
}

eval_results_main <- lapply(hedged_returns_main, evaluate_performance)
eval_results_loo <- lapply(hedged_returns_loo, evaluate_performance)

# 6. Scale-Free Gap Analysis: log(Var_Rolling / Var_Kalman)
get_log_var_ratio <- function(ret_list, valid_c, window = NULL) {
  sapply(valid_c, function(c) {
    dat <- ret_list[[c]]
    if (!is.null(window)) dat <- dat[window]
    v <- apply(dat, 2, var)
    log(v["Rolling"] / v["Kalman"])
  })
}

valid_ccys <- currencies[!sapply(eval_results_main, is.null)]
log_ratio_full <- rbind(Main = get_log_var_ratio(hedged_returns_main, valid_ccys),
                        LOO  = get_log_var_ratio(hedged_returns_loo, valid_ccys))
log_ratio_covid <- rbind(Main = get_log_var_ratio(hedged_returns_main, valid_ccys, "2020-02-01/2020-05-31"),
                         LOO  = get_log_var_ratio(hedged_returns_loo, valid_ccys, "2020-02-01/2020-05-31"))

# 7. Save Metrics
saveRDS(list(
  lb_results = lb_results,
  hedged_returns_main = hedged_returns_main,
  hedged_returns_loo = hedged_returns_loo,
  eval_results_main = eval_results_main,
  eval_results_loo = eval_results_loo,
  log_ratio_full = log_ratio_full,
  log_ratio_covid = log_ratio_covid
), file = "data/04_eval_metrics.rds")

print("Evaluation metrics calculated and saved to data/04_eval_metrics.rds")
