# scripts/05_plots_and_tables.R

# 1. Load Required Packages
library(ggplot2)
library(tidyr)
library(zoo)
library(xts)

# 2. Ensure Output Directories Exist
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# 3. Load Data, Fits, and Metrics
fits <- readRDS("data/03_model_fits.rds")
metrics <- readRDS("data/04_eval_metrics.rds")

# 4. Generate Case Study Plot 1: AUD (Risk-On) COVID-19 Shock
aud_kf <- fits$kf_all[["AUD"]]$oos_hedge$Beta_PC1["2020-01-01/2020-06-30"]
aud_roll <- fits$rolling_oos_all[["AUD"]]$Beta_PC1["2020-01-01/2020-06-30"]

aud_df <- data.frame(Date = index(aud_kf), 
                     Kalman = coredata(aud_kf), 
                     Rolling = coredata(aud_roll))
colnames(aud_df) <- c("Date", "Kalman", "Rolling")

aud_long <- pivot_longer(aud_df, cols = c("Kalman", "Rolling"), names_to = "Model", values_to = "Beta")

p_aud <- ggplot(aud_long, aes(x = Date, y = Beta, color = Model)) +
  geom_line(linewidth = 1) +
  # Use named vectors to strictly enforce color and label mapping
  scale_color_manual(values = c("Kalman" = "#0047AB", 
                                "Rolling" = "#8B0000"),
                     labels = c("Kalman" = "Kalman Filter (State-Space)", 
                                "Rolling" = "252-day Rolling OLS")) +
  labs(title = "AUD PC1 Beta: Kalman vs 252d Rolling OLS (COVID-19)",
       subtitle = "Testing rapid state-space adaptation",
       y = "PC1 Hedge Ratio (Global Dollar)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14))

# Updated to PDF export for Beamer/Paper vector rendering
ggsave("outputs/figures/AUD_COVID_Shock.pdf", plot = p_aud, width = 8, height = 5, bg = "white")

# 5. Generate Case Study Plot 2: JPY (Safe Haven) LOO Anomaly
jpy_kf_main <- fits$kf_all[["JPY"]]$oos_hedge$Beta_PC2["2019-01-01/2025-12-31"]
jpy_kf_loo <- fits$loo_kf_all[["JPY"]]$oos_hedge$Beta_PC2["2019-01-01/2025-12-31"]

jpy_df <- data.frame(Date = index(jpy_kf_main),
                     Main_PCA = coredata(jpy_kf_main),
                     LOO_PCA = coredata(jpy_kf_loo))
colnames(jpy_df) <- c("Date", "Main_PCA", "LOO_PCA")

jpy_long <- pivot_longer(jpy_df, cols = c("Main_PCA", "LOO_PCA"), names_to = "Framework", values_to = "Beta")

p_jpy <- ggplot(jpy_long, aes(x = Date, y = Beta, color = Framework)) +
  geom_line(linewidth = 1) +
  # BUG FIXED: Named vectors explicitly lock the right label to the right data
  scale_color_manual(values = c("Main_PCA" = "#2E8B57", 
                                "LOO_PCA" = "#FF8C00"),
                     labels = c("Main_PCA" = "Standard PCA (Contains JPY)", 
                                "LOO_PCA" = "LOO PCA (Strictly Ex-JPY)")) +
  labs(title = "JPY PC2 Beta",
       subtitle = "Kalman Filter's dynamic estimation with and without LOO",
       y = "PC2 Hedge Ratio (Carry)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14))

# Updated to PDF export for Beamer/Paper vector rendering
ggsave("outputs/figures/JPY_LOO_Anomaly.pdf", plot = p_jpy, width = 8, height = 5, bg = "white")

# 6. Extract and Export Master Tables to CSV
currencies <- names(metrics$eval_results_main)

# Function to extract metrics from Main PCA results
extract_main_metric <- function(metric_name) {
  t(sapply(currencies, function(ccy) metrics$eval_results_main[[ccy]][metric_name, ]))
}

# Function to extract metrics from LOO PCA results
extract_loo_metric <- function(metric_name) {
  t(sapply(currencies, function(ccy) metrics$eval_results_loo[[ccy]][metric_name, ]))
}

# --- Export Main (Standard PCA) Metrics ---
write.csv(extract_main_metric("Variance_Reduction_Full"), 
          "outputs/tables/Main_Variance_Reduction_Full.csv")
write.csv(extract_main_metric("Variance_Reduction_COVID"), 
          "outputs/tables/Main_Variance_Reduction_COVID.csv")
write.csv(extract_main_metric("VaR_99_Full"), 
          "outputs/tables/Main_VaR_99_Full.csv")
write.csv(extract_main_metric("Max_Drawdown_COVID"), 
          "outputs/tables/Main_Max_Drawdown_COVID.csv")

# --- Export LOO (Leave-One-Out) Metrics ---
write.csv(extract_loo_metric("Variance_Reduction_Full"), 
          "outputs/tables/LOO_Variance_Reduction_Full.csv")
write.csv(extract_loo_metric("Variance_Reduction_COVID"), 
          "outputs/tables/LOO_Variance_Reduction_COVID.csv")
write.csv(extract_loo_metric("VaR_99_Full"), 
          "outputs/tables/LOO_VaR_99_Full.csv")
write.csv(extract_loo_metric("Max_Drawdown_COVID"), 
          "outputs/tables/LOO_Max_Drawdown_COVID.csv")

# Export Ljung-Box Diagnostics
write.csv(round(metrics$lb_results, 4), 
          "outputs/tables/Ljung_Box_Diagnostics.csv")

print("Plots generated and saved to outputs/figures/")
print("Main and LOO Master tables generated and saved to outputs/tables/")