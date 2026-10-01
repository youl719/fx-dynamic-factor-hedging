# 1. Load Required Packages
# install.packages(c("quantmod", "dplyr", "tidyr", "zoo"))
library(quantmod)
library(dplyr)
library(tidyr)
library(zoo)

# 2. Define FRED Series Tickers for G10 Currencies
# Note: We exclude USD as it's the numeraire. We pull 9 pairs.
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

# 4. Clean Missing Data (US Holidays)
# Carry forward the previous day's closing price
fx_prices <- na.locf(fx_prices) 

# 5. Standardize Quoting & Calculate Log Returns
# Invert the required pairs so all represent 1 Unit of Foreign Currency = X USD
fx_prices$JPY <- 1 / fx_prices$JPY
fx_prices$CAD <- 1 / fx_prices$CAD
fx_prices$CHF <- 1 / fx_prices$CHF
fx_prices$SEK <- 1 / fx_prices$SEK
fx_prices$NOK <- 1 / fx_prices$NOK

# Calculate daily log returns
fx_returns <- diff(log(fx_prices))
fx_returns <- na.omit(fx_returns) # Drop the initial NA row

# 6. Split into In-Sample (Train) and Out-of-Sample (Test)
train_window <- paste0("2005-01-01/", "2018-12-31")
test_window <- paste0("2019-01-01/", Sys.Date())

train_returns <- fx_returns[train_window]
test_returns <- fx_returns[test_window]

# 7. Pure PCA on the Training Data
# scale. = TRUE standardizes the variances (uses correlation matrix instead of covariance)
pca_train <- prcomp(train_returns, center = TRUE, scale. = TRUE)

# Optional: View variance explained. PC1 and PC2 usually explain > 70% of G10 variance.
summary(pca_train)

# Extract Training Factors (In-Sample)
factors_train <- pca_train$x[, 1:2]
colnames(factors_train) <- c("PC1_Dollar", "PC2_Carry")
factors_train <- as.xts(factors_train, order.by = index(train_returns))

# 8. Out-of-Sample Projection (Zero Look-Ahead Bias)
# Scale the test data using ONLY the training mean and standard deviation
train_mean <- pca_train$center
train_sd <- pca_train$scale

test_returns_scaled <- scale(test_returns, center = train_mean, scale = train_sd)
test_returns_scaled[is.na(test_returns_scaled)] <- 0 # Failsafe for any stray NAs

# Multiply scaled test returns by the training eigenvectors (rotation matrix)
factors_test <- test_returns_scaled %*% pca_train$rotation[, 1:2]
colnames(factors_test) <- c("PC1_Dollar", "PC2_Carry")

# Convert back to an xts object for easy modeling
factors_test <- as.xts(factors_test, order.by = index(test_returns))

# Quick plot of the OOS COVID-19 shock on PC1
plot(factors_test["2020-01-01/2020-06-30"]$PC1_Dollar, 
     main = "OOS PC1 (Dollar Factor) During COVID Shock", 
     col = "blue")

# 1. Combine Training Returns and Factors
# Assuming train_returns and factors_train are already xts objects in your environment
train_data <- merge(train_returns, factors_train)
train_data <- na.omit(train_data) # Ensure no missing data at the very beginning

# 2. Static OLS: Full Sample (2005-2018)
# Loop through all 9 currencies to generate a master table of static betas
currencies <- c("EUR", "GBP", "AUD", "NZD", "JPY", "CAD", "CHF", "SEK", "NOK")

static_betas <- sapply(currencies, function(ccy) {
  # Construct the formula dynamically (e.g., "AUD ~ PC1_Dollar + PC2_Carry")
  formula_str <- paste(ccy, "~ PC1_Dollar + PC2_Carry")
  model <- lm(as.formula(formula_str), data = train_data)
  return(coef(model))
})

# Transpose for readability
static_betas_table <- t(static_betas)
print("Static OLS Betas (2005-2018):")
print(static_betas_table)

# 3. Rolling-Window OLS: 252 Days (Example: AUD)
# Define a custom function to pass into rollapply
run_rolling_ols <- function(data_subset) {
  # data_subset[, 1] is the target currency, [, 2] is PC1, [, 3] is PC2
  model <- lm(data_subset[, 1] ~ data_subset[, 2] + data_subset[, 3])
  return(coef(model))
}

# Isolate the AUD dataset
aud_rolling_data <- train_data[, c("AUD", "PC1_Dollar", "PC2_Carry")]

# Execute rolling regression
# width = 252 (approx. 1 trading year)
# by.column = FALSE ensures the function takes the whole matrix block at once
rolling_betas_aud <- rollapply(
  data = aud_rolling_data,
  width = 252,
  FUN = run_rolling_ols,
  by.column = FALSE,
  align = "right" 
)

# Clean up the output and convert back to an xts time-series object
colnames(rolling_betas_aud) <- c("Alpha", "Beta_PC1", "Beta_PC2")
rolling_betas_aud <- as.xts(rolling_betas_aud)

# 4. Visualizing the Instability (Fix for xts plotting quirks)
# Convert the rolling beta to a standard numeric vector for a foolproof plot
plot(as.numeric(rolling_betas_aud$Beta_PC1), 
     type = "l", 
     col = "darkred", 
     main = "AUD Rolling Beta vs. Static Beta Assumption",
     ylab = "PC1 Beta", 
     xlab = "Trading Days (2005-2018)")

# The blue line should now perfectly overlay across the chart
abline(h = static_betas_table["AUD", "PC1_Dollar"], col = "blue", lwd = 2, lty = 2)