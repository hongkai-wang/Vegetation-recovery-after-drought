library(ranger)
library(dplyr)

predictors <- c("drought_duration", "drought_intensity",
                "temp_mean", "pre_mean", "pet_mean", "sm_mean",
                "ssrd_mean", "kndvi_mean")

model_recovery_rate <- results_combined_final %>%
  select(lai_recovery_rate, wue_recovery_rate, all_of(predictors)) %>%
  na.omit() 

set.seed(123)
sample_n <- min(1000, nrow(model_recovery_rate))  # 最多1000条
model_sample <- sample_n(model_recovery_rate, sample_n)

# 模型A：LAI恢复速率
model_lai_rate <- ranger(
  formula = lai_recovery_rate ~ .,
  data = model_sample %>% select(lai_recovery_rate, all_of(predictors)),
  importance = "impurity",      
  num.trees = 200,             
  write.forest = TRUE
)

# 模型B：WUE恢复速率
model_wue_rate <- ranger(
  formula = wue_recovery_rate ~ .,
  data = model_sample %>% select(wue_recovery_rate, all_of(predictors)),
  importance = "impurity",
  num.trees = 200,
  write.forest = TRUE
)

# === 6. 输出变量重要性 ===
cat("\n🟢 LAI Recovery Rate - Variable Importance:\n")
print(sort(model_lai_rate$variable.importance, decreasing = TRUE))

cat("\n🔵 WUE Recovery Rate - Variable Importance:\n")
print(sort(model_wue_rate$variable.importance, decreasing = TRUE))

pred_lai <- predict(model_lai_rate,
                    data = model_sample %>% select(all_of(predictors)))$predictions

pred_wue <- predict(model_wue_rate,
                    data = model_sample %>% select(all_of(predictors)))$predictions

true_lai <- model_sample$lai_recovery_rate
true_wue <- model_sample$wue_recovery_rate

# R² 和 RMSE 计算函数
r2 <- function(true, pred) 1 - sum((true - pred)^2) / sum((true - mean(true))^2)
rmse <- function(true, pred) sqrt(mean((true - pred)^2))

# 输出
cat("🟢 LAI Recovery Rate:\nR² =", round(r2(true_lai, pred_lai), 3),
    "  RMSE =", round(rmse(true_lai, pred_lai), 3), "\n")

cat("🔵 WUE Recovery Rate:\nR² =", round(r2(true_wue, pred_wue), 3),
    "  RMSE =", round(rmse(true_wue, pred_wue), 3), "\n")


library(ggplot2)

# 整理变量重要性数据
importance_df_rate <- bind_rows(
  data.frame(Variable = names(model_lai_rate$variable.importance),
             Importance = as.numeric(model_lai_rate$variable.importance),
             Target = "LAI Recovery Rate"),
  data.frame(Variable = names(model_wue_rate$variable.importance),
             Importance = as.numeric(model_wue_rate$variable.importance),
             Target = "WUE Recovery Rate")
)

# 绘图
ggplot(importance_df_rate, aes(x = reorder(Variable, Importance), y = Importance, fill = Target)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Variable Importance (Random Forest)",
       x = "Predictor", y = "Importance") +
  theme_minimal()



#建模 LAI 和 WUE 滞后期lag
predictors <- c("drought_duration", "drought_intensity",
                "temp_mean", "pre_mean", "pet_mean", "sm_mean",
                "ssrd_mean", "kndvi_mean")

model_lag_data <- results_combined_final %>%
  select(lai_lag_to_positive, wue_lag_to_positive, all_of(predictors)) %>%
  na.omit()  

set.seed(123)
sample_n <- min(1000, nrow(model_lag_data))  # 最多1000条
model_lag_sample <- sample_n(model_lag_data, sample_n)

#LAI恢复滞后时间
model_lai_lag <- ranger(
  formula = lai_lag_to_positive ~ .,
  data = model_lag_sample %>% select(lai_lag_to_positive, all_of(predictors)),
  importance = "impurity",    
  num.trees = 200,             
  write.forest = TRUE
)

#模型B：WUE恢复滞后时间 
model_wue_lag <- ranger(
  formula = wue_lag_to_positive ~ .,
  data = model_lag_sample %>% select(wue_lag_to_positive, all_of(predictors)),
  importance = "impurity",
  num.trees = 200,
  write.forest = TRUE
)

# === 5. 输出变量重要性 ===
cat("\n🟢 LAI Lag to Positive - Variable Importance:\n")
print(sort(model_lai_lag$variable.importance, decreasing = TRUE))

cat("\n🔵 WUE Lag to Positive - Variable Importance:\n")
print(sort(model_wue_lag$variable.importance, decreasing = TRUE))

pred_lai_lag <- predict(model_lai_lag,
                        data = model_lag_sample %>% select(all_of(predictors)))$predictions

pred_wue_lag <- predict(model_wue_lag,
                        data = model_lag_sample %>% select(all_of(predictors)))$predictions

true_lai_lag <- model_lag_sample$lai_lag_to_positive
true_wue_lag <- model_lag_sample$wue_lag_to_positive

# 输出评估指标
cat("🟢 LAI Lag to Positive:\nR² =", round(r2(true_lai_lag, pred_lai_lag), 3),
    "  RMSE =", round(rmse(true_lai_lag, pred_lai_lag), 3), "\n")

cat("🔵 WUE Lag to Positive:\nR² =", round(r2(true_wue_lag, pred_wue_lag), 3),
    "  RMSE =", round(rmse(true_wue_lag, pred_wue_lag), 3), "\n")

# === 7. 可视化变量重要性 ===
library(ggplot2)

# 整理变量重要性数据
importance_df_lag <- bind_rows(
  data.frame(Variable = names(model_lai_lag$variable.importance),
             Importance = as.numeric(model_lai_lag$variable.importance),
             Target = "LAI Lag to Positive"),
  data.frame(Variable = names(model_wue_lag$variable.importance),
             Importance = as.numeric(model_wue_lag$variable.importance),
             Target = "WUE Lag to Positive")
)

# 绘图
ggplot(importance_df_lag, aes(x = reorder(Variable, Importance), y = Importance, fill = Target)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Variable Importance for Lag to Positive (Random Forest)",
       x = "Predictor", y = "Importance") +
  theme_minimal() +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62"))




