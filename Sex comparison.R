library(dplyr)
library(ggplot2)

# Relabel surveys (PPS3 -> 2024, PPS4 -> 2025); 2024 will sit on top
sex_df <- sex_df %>%
  mutate(
    survey = recode(as.character(survey), "PPS3" = "2024", "PPS4" = "2025"),
    survey = factor(survey, levels = c("2025", "2024"))
  )

p <- ggplot(sex_df, aes(x = pct, y = survey, fill = sex)) +
  geom_col(width = 0.5, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 4.2) +
  scale_fill_manual(values = c(Female = "#6250d6", Male = "#1baf7a"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(add = c(0, 4))          # room on the right so 100% isn't clipped
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "Sex split of patients on antimicrobials",
    subtitle = "2024 (n = 2,102) vs 2025 (n = 1,985)",
    x = NULL, y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor     = element_blank(),
    panel.grid.major.y   = element_blank(),
    panel.grid.major.x   = element_line(colour = "#E8E8E8"),
    axis.text.x          = element_text(size = 12, colour = "#333333",
                                        margin = margin(t = 6)),
    axis.text.y          = element_text(size = 13, face = "bold", colour = "#333333"),
    axis.line.x          = element_line(colour = "#999999"),
    axis.ticks.x         = element_line(colour = "#999999"),
    axis.ticks.length.x  = unit(0.2, "cm"),
    plot.title           = element_text(face = "bold"),
    legend.position      = "top",
    legend.justification = "left",
    plot.margin          = margin(15, 25, 10, 15)
  )

print(p)

# 5. Save -------------------------------------------------------------------
ggsave("02pps_sex_split.png", p, width = 7, height = 3, dpi = 300)