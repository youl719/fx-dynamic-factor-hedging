# scripts/01_data_and_factors.R

# 1. Load Required Packages
library(quantmod)
library(zoo)

# 2. Define FRED Series Tickers for G10 Currencies
tickers <- c(
  EUR = "DEXUSEU", GBP = "DEXUSUK", AUD = "DEXUSAL", NZD = "DEXUSNZ",
  JPY = "DEXJPUS", CAD = "DEXCAUS", CHF = "DEXSZUS", SEK = "DEXSDUS", NOK = "DEXNOUS"
)

# 3. Download Data via quantmod
getSymbols(tickers, src = "FRED", env = globalenv())
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

fx_returns <- na.omit(diff(log(fx_prices)))

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

scale_fx <- function(x) scale(x, center = pca_train$center, scale = pca_train$scale)
factors_train <- as.xts(scale_fx(train_returns) %*% rot, order.by = index(train_returns))
colnames(factors_train) <- c("PC1_Dollar", "PC2_Carry")

# 8. Out-of-Sample Projection
test_returns_scaled <- scale_fx(test_returns)
factors_test <- as.xts(test_returns_scaled %*% rot, order.by = index(test_returns))
colnames(factors_test) <- c("PC1_Dollar", "PC2_Carry")

# 9. Combine Train and Test for Continuous Modeling
all_factors <- rbind(factors_train, factors_test)
all_returns <- rbind(train_returns, test_returns)
all_data <- na.omit(merge(all_returns, all_factors))

# NEW: Save the finalized environment for the next scripts
saveRDS(list(
  all_data = all_data,
  train_returns = train_returns,
  test_returns = test_returns,
  train_window = train_window,
  test_window = test_window,
  currencies = names(tickers)
), file = "data/01_cleaned_data.rds")

print("Data cleaning and main PCA complete. Objects saved to data/01_cleaned_data.rds")
