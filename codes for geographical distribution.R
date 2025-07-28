library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# === 读取世界地图数据 ===
world_land <- ne_countries(scale = "medium", returnclass = "sf")

# === 设置类型与配色 ===
results_processed <- results_processed %>%
  mutate(recovery_type_label = factor(
    recovery_type_label,
    levels = c("Type1 Full Recovery", 
               "Type2 Structure Dominant", 
               "Type3 Function Dominant", 
               "Type4 Failed Recovery")
  ))

custom_colors <- c(
  "Type1 Full Recovery" = "#1f78b4",       # 蓝色
  "Type2 Structure Dominant" = "#33a02c",  # 绿色
  "Type3 Function Dominant" = "#ff7f00",   # 橙色
  "Type4 Failed Recovery" = "#e31a1c"      # 红色
)

# === 绘图 ===
combined_plot <- ggplot() +
  geom_point(data = results_processed,
             aes(x = lon, y = lat, color = recovery_type_label),
             size = 1.2, alpha = 0.8, shape = 16) +
  geom_sf(data = world_land, fill = NA, color = "black", linewidth = 0.3) +
  
  scale_color_manual(
    values = custom_colors,
    name = NULL,
    guide = guide_legend(
      direction = "horizontal",
      nrow = 1,
      label.position = "bottom",
      title.position = "top",
      override.aes = list(size = 5),    
      label.hjust = 0.5,                
      label.vjust = 1.2,                
      label.theme = element_text(margin = margin(t = 6, r = 10, b = 0, l = 10)) 
    )
  ) +
  
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 90), expand = FALSE) +
  
  labs(
    title = "",
    x = "Longitude", y = "Latitude"
  ) +
  
  theme_void(base_size = 16, base_family = "Times New Roman") +
  
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 3),
    axis.text = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold", angle = 90),
    
    # 边框和网格
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 14) 
  )

# === 保存图像 ===
ggsave("Fig_Global_Recovery_Types_Map_LegendBottom_Enhanced.png",
       combined_plot, width = 14.5, height = 7, dpi = 600, bg = "white")

# === 显示图像 ===
print(combined_plot)


