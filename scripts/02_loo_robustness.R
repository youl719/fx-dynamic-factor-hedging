# scripts/02_loo_robustness.R

# 1. Load Required Packages
library(zoo)

# 2. Load the cleaned data from Step 1
env_data <- readRDS("data/01_cleaned_data.rds")
all_data <- env_data$all_data
train_returns <- env_data$train_returns
test_returns <- env_data$test_returns
train_window <- env_data$train_window
test_window <- env_data$test_window
currencies <- env_data$currencies

# 3. Leave-One-Out (LOO) PCA, Data Structuring, & Correlation Alignment
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
  
  cors <- c(PC1 = cor(as.numeric(f[train_window][, 1]), as.numeric(main_tr[, 1])),
            PC2 = cor(as.numeric(f[train_window][, 2]), as.numeric(main_tr[, 2])))
  
  f[, 1] <- sign(cors[1]) * f[, 1]
  f[, 2] <- sign(cors[2]) * f[, 2]
  
  dat <- na.omit(merge(rbind(train_returns[, ccy], test_returns[, ccy]), f))
  colnames(dat) <- c(ccy, "PC1_Dollar", "PC2_Carry")
  
  mod <- lm(as.formula(paste(ccy, "~ PC1_Dollar + PC2_Carry")),
            data = as.data.frame(dat[train_window]))
  
  list(data = dat, static = coef(mod), cors = abs(cors), eigenvectors = rot)
}

# 4. Execute LOO for all currencies
loo_objs <- setNames(lapply(currencies, build_loo), currencies)
loo_data_list <- lapply(loo_objs, `[[`, "data")
loo_static_betas_table <- t(sapply(loo_objs, `[[`, "static"))
loo_cors <- round(t(sapply(loo_objs, `[[`, "cors")), 3)
loo_eigenvectors <- lapply(loo_objs, `[[`, "eigenvectors")

# 5. Save LOO outputs for the modeling script
saveRDS(list(
  loo_data_list = loo_data_list,
  loo_static_betas_table = loo_static_betas_table,
  loo_cors = loo_cors,
  loo_eigenvectors = loo_eigenvectors
), file = "data/02_loo_data.rds")
