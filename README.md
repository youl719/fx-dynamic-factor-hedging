# Dynamic vs. Static Factor Hedging in FX: A Quantitative Analysis of G10 Currencies

An empirical econometrics study comparing dynamic state-space factor models (estimated via the Kalman filter) against static factor benchmarks for hedging currency risk across G10 FX portfolios.

## Overview & Motivation
Traditional static factor hedging relies on constant factor loadings estimated over historical windows, which break down during macroeconomic regime shifts and volatility spikes (such as the COVID-19 market shock). This project evaluates whether time-varying factor exposures modeled via a state-space Kalman filter improve out-of-sample variance reduction, tail risk (99% VaR), and maximum drawdown compared to static factor models, alongside leave-one-out (LOO) cross-currency robustness checks.

## Methodology
* **Data:** Daily/weekly G10 currency returns and macroeconomic/market risk factors (see `data/`).
* **Static Benchmark:** Ordinary Least Squares (OLS) / rolling-window factor regressions estimating constant hedge ratios.
* **Dynamic Model:** Rolling window OLS and Linear Gaussian state-space model where latent factor loadings follow a random walk process, estimated via Maximum Likelihood and the **Kalman Filter**.
* **Evaluation & Stress Testing:** Out-of-sample portfolio variance reduction, 99% Value-at-Risk (VaR), and Maximum Drawdown across both the full sample and the COVID-19 stress period.
## Key Findings
* **Dynamic vs Static factor hedging:** Outperformance of dynamic models to the static baseline was consistent across all risk metrics
* **Comparison between dynamic models:** Kalman filter performed better during COVID crisis on the safe haven assets but underperformed in the remaining riskier currencies.

> **Full Paper:** See [`outputs/Financial_econometrics_report.pdf`](outputs/Financial_econometrics_report.pdf) 

## Repository Structure
```text
fx-dynamic-factor-hedging/
├── data/                    # Processed datasets and intermediate .rds evaluation objects
├── scripts/                 # Modular R pipeline
│   ├── 01_data_and_factors.R              # Data cleaning and factor construction
│   ├── 02_loo_robustness.R  # Leave-one-out currency sensitivity analysis
│   ├── 04_evaluation_metrics.R # Out-of-sample variance, 99% VaR, and drawdown metrics
│   └── 05_plots_and_tables.R   # Generation of publication figures and summary tables
├── outputs/
│   ├── Financial_econometrics_report.pdf  # Final research report
│   ├── figures/             # Crisis and LOO diagnostic plots (AUD, JPY)
  │   └── tables/              # Main and LOO performance comparison tables (.csv)
└── .gitignore

