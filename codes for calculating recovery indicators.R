library(ncdf4)
library(terra)
library(tibble)
library(dplyr)
library(lubridate)
# === 加载 LAI Z-score NetCDF 数据 ===
nc_lai <- nc_open("E:/data_paper2/Global/LAI_global_050res_zscore_200003_202212.nc")
lai_array <- ncvar_get(nc_lai, "LAI_zscore")  
lai_array <- aperm(lai_array, c(2, 1, 3))  
# === 加载 WUE Z-score NetCDF 数据 ===
nc_wue <- nc_open("E:/data_paper2/Global/WUE_global_050res_zscore_200003_202212.nc")
wue_array <- ncvar_get(nc_wue, "WUE_zscore") 
nc_close(nc_wue)
wue_array <- aperm(wue_array, c(2, 1, 3)) 

# === 时间序列定义 ===
time_seq <- seq(as.Date("2000-03-01"), by = "month", length.out = dim(lai_array)[3])
# === 初始化结果列表 ===
results_list <- list()
# === 提取每个事件的恢复期 LAI 和 WUE ===
for (i in seq_len(nrow(SPEI_event_df))) {
  row <- SPEI_event_df[i, ]
  lat_idx <- row$lat_idx
  lon_idx <- row$lon_idx
  recovery_start <- which(time_seq == as.Date(row$recovery_start_date))
  recovery_end <- which(time_seq == as.Date(row$recovery_end_date))
  
  if (length(recovery_start) == 0 || length(recovery_end) == 0) next
  # 提取 LAI 序列
  lai_ts <- lai_array[lat_idx, lon_idx, recovery_start:recovery_end]
  if (all(is.na(lai_ts))) next
  # 提取 WUE 序列
  wue_ts <- wue_array[lat_idx, lon_idx, recovery_start:recovery_end]
  # 计算恢复时长与类型
  recovery_duration <- recovery_end - recovery_start + 1
  recovery_type <- case_when(
    recovery_duration <= 6 ~ "fast",
    recovery_duration <= 12 ~ "mid",
    TRUE ~ "long"
  )
  result_row <- tibble(
    event_id = i,
    lat = row$lat,
    lon = row$lon,
    lat_idx = lat_idx,
    lon_idx = lon_idx,
    drought_start_date = row$drought_start_date,
    drought_end_date = row$drought_end_date,
    recovery_start_date = row$recovery_start_date,
    recovery_end_date = row$recovery_end_date,
    drought_duration = row$drought_duration,
    drought_intensity = row$drought_intensity,
    recovery_duration = recovery_duration,
    recovery_type = recovery_type,
    recovery_dates = list(time_seq[recovery_start:recovery_end]),
    lai_zscore = list(as.numeric(lai_ts)),
    wue_zscore = list(as.numeric(wue_ts))
  )
  
  results_list[[i]] <- result_row
}

# 合并所有结果为数据框 
results_LAIWUE <- bind_rows(results_list)


library(dplyr)
library(purrr)
library(tibble)

is_increasing <- function(ts) {
  if (all(is.na(ts)) || length(na.omit(ts)) < 2) return(FALSE)
  time_index <- seq_along(ts)
  fit <- lm(ts ~ time_index)
  return(coef(fit)[2] > 0)
}

results_LAIWUE_increased <- results_LAIWUE %>%
  rowwise() %>%
  filter(
    !is.null(lai_zscore) && !is.null(wue_zscore),
    length(lai_zscore) > 0,
    length(wue_zscore) > 0,
    any(lai_zscore < 0, na.rm = TRUE),  # 有低谷即判定为受干扰
    any(wue_zscore < 0, na.rm = TRUE),
    is_increasing(lai_zscore),
    is_increasing(wue_zscore)
  ) %>%
  ungroup()


#恢复指标的计算和恢复类型的分类
library(dplyr)
library(purrr)

# 判断首次出现“连续 N 月 Z-score ≥ 阈值”的位置
get_lag_to_positive_strict <- function(ts, min_consecutive = 4, threshold = 0) {
  n <- length(ts)
  if (n < min_consecutive || all(is.na(ts))) return(NA)
  
  for (i in seq_len(n - min_consecutive + 1)) {
    window <- ts[i:(i + min_consecutive - 1)]
    if (all(window >= threshold, na.rm = TRUE)) {
      return(i)
    }
  }
  return(NA)
}

# 恢复指标提取函数：整合 stricter 滞后期判定
extract_metrics <- function(zs, min_consecutive = 4, threshold = 0) {
  n <- length(zs)
  if (n == 0 || all(is.na(zs))) {
    return(list(min = NA, max = NA, recovery_rate = NA, lag_to_positive = NA,
                early_mean = NA, mid_mean = NA, late_mean = NA))
  }
  
  min_val <- min(zs, na.rm = TRUE)
  min_idx <- which.min(zs)
  
  if (min_idx < n) {
    post_zs <- zs[(min_idx + 1):n]
    if (all(is.na(post_zs))) {
      rate <- NA
    } else {
      post_max_val <- max(post_zs, na.rm = TRUE)
      post_max_rel_idx <- which.max(post_zs)
      rate <- (post_max_val - min_val) / post_max_rel_idx
    }
    lag <- get_lag_to_positive_strict(post_zs, min_consecutive = min_consecutive, threshold = threshold)
  } else {
    rate <- NA
    lag <- NA
  }
  
  early_mean <- mean(zs[1:min(6, n)], na.rm = TRUE)
  mid_mean   <- mean(zs[7:min(12, n)], na.rm = TRUE)
  late_mean  <- mean(zs[13:min(24, n)], na.rm = TRUE)
  
  return(list(min = min_val, max = max(zs, na.rm = TRUE), recovery_rate = rate,
              lag_to_positive = lag, early_mean = early_mean,
              mid_mean = mid_mean, late_mean = late_mean))
}


#分类函数（基于转正滞后期）
get_recovery_type <- function(lai_zs, wue_zs) {
  lai_lag <- extract_metrics(lai_zs)$lag_to_positive
  wue_lag <- extract_metrics(wue_zs)$lag_to_positive
  
  lai_pos <- !is.na(lai_lag)
  wue_pos <- !is.na(wue_lag)
  
  if (lai_pos && wue_pos) {
    return("Type1 Full Recovery")
  } else if (lai_pos && !wue_pos) {
    return("Type2 Structure Dominant")
  } else if (!lai_pos && wue_pos) {
    return("Type3 Function Dominant")
  } else {
    return("Type4 Failed Recovery")
  }
}

# 应用于主数据框
results_processed <- results_LAIWUE_increased %>%
  rowwise() %>%
  mutate(
    lai_metrics = list(extract_metrics(lai_zscore)),
    lai_min = lai_metrics$min,
    lai_max = lai_metrics$max,
    lai_recovery_rate = lai_metrics$recovery_rate,
    lai_lag_to_positive = lai_metrics$lag_to_positive,
    lai_early_mean = lai_metrics$early_mean,
    lai_mid_mean = lai_metrics$mid_mean,
    lai_late_mean = lai_metrics$late_mean,
    
    wue_metrics = list(extract_metrics(wue_zscore)),
    wue_min = wue_metrics$min,
    wue_max = wue_metrics$max,
    wue_recovery_rate = wue_metrics$recovery_rate,
    wue_lag_to_positive = wue_metrics$lag_to_positive,
    wue_early_mean = wue_metrics$early_mean,
    wue_mid_mean = wue_metrics$mid_mean,
    wue_late_mean = wue_metrics$late_mean,
    
    recovery_type_label = get_recovery_type(lai_zscore, wue_zscore)
  ) %>%
  ungroup() %>%
  dplyr::select(-lai_metrics, -wue_metrics)

# 可视化或统计输出
table(results_processed$recovery_type_label)


#绘制数量占比
library(ggplot2)
library(dplyr)

# === 1. 计算占比 ===
results_df <- as.data.frame(st_drop_geometry(results_processed))

# 统计恢复类型的数量和占比
type_counts <- results_df %>%
  group_by(recovery_type_label) %>%
  summarise(n = n()) %>%
  mutate(Percent = n / sum(n) * 100)
# 2. 设置颜色
custom_colors <- c(
  "Type1 Full Recovery" = "#3B528B", 
  "Type2 Structure Dominant" = "#4393C3",  
  "Type3 Function Dominant" = "#7BCCC4",  
  "Type4 Failed Recovery" = "#FDE725" 
)

p1 <- ggplot(type_counts, aes(x = reorder(recovery_type_label, -Percent), 
                             y = Percent, fill = recovery_type_label)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = custom_colors, name = "Recovery Type") +
  labs(
    title = "Proportions of Vegetation Drought Recovery Types",
    x = NULL, y = "Percentage (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 13, angle = 15, hjust = 1),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",  # 图例在图中不需要
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Fig_Recovery_Type_Proportions_Bar.png",
       plot = p, width = 8, height = 5, dpi = 300, bg = "white")

print(p1)

#不同恢复类型的恢复指标箱线图
names(results_processed)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggsci)
library(ggpubr)
# === 1. 添加 Range 指标到原始数据 ===
results_processed <- results_processed %>%
  mutate(
    lai_range = lai_max - lai_min,
    wue_range = wue_max - wue_min
  )
# === 2. 数据整理（加入 Range）===
plot_data <- results_processed %>%
  select(recovery_type_label,
         lai_range, wue_range,
         lai_recovery_rate, wue_recovery_rate,
         lai_lag_to_positive, wue_lag_to_positive,
         lai_early_mean, wue_early_mean,
         lai_mid_mean, wue_mid_mean,
         lai_late_mean, wue_late_mean) %>%
  pivot_longer(-recovery_type_label, names_to = "variable", values_to = "value") %>%
  mutate(
    recovery_type_label = recode(recovery_type_label,
                                 "Type1 Full Recovery" = "Type1",
                                 "Type2 Structure Dominant" = "Type2",
                                 "Type3 Function Dominant" = "Type3",
                                 "Type4 Failed Recovery" = "Type4"),
    group = ifelse(grepl("^lai", variable), "LAI", "WUE"),
    indicator = case_when(
      grepl("range", variable) ~ "Range",
      grepl("recovery_rate", variable) ~ "Recovery Rate",
      grepl("lag_to_positive", variable) ~ "Lag to Positive",
      grepl("early_mean", variable) ~ "Early Mean",
      grepl("mid_mean", variable) ~ "Mid Mean",
      grepl("late_mean", variable) ~ "Late Mean"
    ),
    group = factor(group, levels = c("LAI", "WUE")),
    indicator = factor(indicator, levels = c(
      "Range", "Recovery Rate", "Lag to Positive",
      "Early Mean", "Mid Mean", "Late Mean"
    )),
    recovery_type_label = fct_relevel(recovery_type_label, "Type1", "Type2", "Type3", "Type4")
  )

#3. 显著性检验
valid_combos <- plot_data %>%
  group_by(indicator, recovery_type_label, group) %>%
  summarise(n = sum(!is.na(value)), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  filter(LAI >= 2 & WUE >= 2) %>%
  select(indicator, recovery_type_label)

plot_data_valid <- plot_data %>%
  semi_join(valid_combos, by = c("indicator", "recovery_type_label"))

stat_tests <- compare_means(
  value ~ group,
  group.by = c("indicator", "recovery_type_label"),
  data = plot_data_valid,
  method = "t.test"
)

stat_tests <- stat_tests %>%
  left_join(
    plot_data_valid %>%
      group_by(indicator, recovery_type_label) %>%
      summarise(y.position = max(value, na.rm = TRUE) * 1.05, .groups = "drop"),
    by = c("indicator", "recovery_type_label")
  )

# 4. 自定义 Times New Roman 主题
theme_sci <- function(base_size = 13) {
  theme_minimal(base_size = base_size, base_family = "Times New Roman") +
    theme(
      text = element_text(family = "Times New Roman", color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(size = rel(0.9), color = "black"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey60", fill = NA),
      strip.background = element_rect(fill = "#f0f0f0", color = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(1.5), family = "Times New Roman"), 
      plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.3))
    )
}
# 5. 绘图
p2 <- ggplot(plot_data, aes(x = recovery_type_label, y = value, fill = group)) +
  geom_violin(
    aes(group = interaction(recovery_type_label, group)),
    trim = FALSE, scale = "width", alpha = 0.5, color = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_boxplot(
    aes(group = interaction(recovery_type_label, group)),
    width = 0.15, outlier.shape = NA, alpha = 0.9, color = "black",
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~indicator, scales = "free_y", nrow = 2, ncol = 3) +  
  scale_fill_manual(values = c("LAI" = "#2CA02C", "WUE" = "#1F77B4")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 2)) +  
  geom_text(
    data = stat_tests,
    aes(x = recovery_type_label, y = y.position * 1.1, label = p.signif),
    inherit.aes = FALSE,
    size = 4.2, fontface = "bold", color = "black", family = "Times New Roman"
  ) +
  labs(
    x = "Recovery Type",
    y = NULL
  ) +
  theme_sci()
ggsave(
  filename = "LAI_WUE_Recovery_Indicators_with_Range.png",
  plot = p2,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)
