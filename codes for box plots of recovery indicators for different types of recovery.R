library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggsci)
library(ggpubr)
results_processed <- results_processed %>%
  mutate(
    lai_range = lai_max - lai_min,
    wue_range = wue_max - wue_min
  )

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

#显著性检验
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

#绘图
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


print(p2)

ggsave(
  filename = "LAI_WUE_Recovery_Indicators_with_Range.png",
  plot = p2,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)
