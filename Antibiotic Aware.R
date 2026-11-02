# ============================================================
# 8. E. coli resistance trends and MAR index
# ============================================================

library(tidyverse)
library(scales)
library(patchwork)

# 1. Resistance by antibiotic: horizontal grouped bar plot
p8a <- ggplot(
  trend,
  aes(
    x = pct_R,
    y = reorder(antibiotic, pct_R),
    fill = round
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_R)),
    position = position_dodge(width = 0.75),
    hjust = -0.15,
    size = 3.2
  ) +
  scale_x_continuous(
    labels = label_percent(scale = 1),
    limits = c(0, 110),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = c("Round 1" = "#2878B5", "Round 2" = "#E07A3F")
  ) +
  labs(
    title = "Antimicrobial resistance in E. coli",
    subtitle = "Resistance prevalence by antibiotic and sampling round",
    x = "Resistant isolates (%)",
    y = NULL,
    fill = "Sampling round"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey35"),
    axis.text.y = element_text(color = "grey20"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.margin = margin(10, 25, 10, 10)
  )

# Save resistance plot
save_plot(p8a, "08_ecoli_resistance_trend", w = 10, h = 7)


# 2. Mean MAR index by sampling round
mar_summary <- mar %>%
  group_by(round) %>%
  summarise(
    n_isolates = n(),
    mean_MAR = mean(mar, na.rm = TRUE),
    median_MAR = median(mar, na.rm = TRUE),
    .groups = "drop"
  )

p8b <- ggplot(
  mar_summary,
  aes(x = mean_MAR, y = reorder(round, mean_MAR), fill = round)
) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.3f", mean_MAR)),
    hjust = -0.2,
    size = 4
  ) +
  scale_x_continuous(
    limits = c(0, max(mar_summary$mean_MAR, na.rm = TRUE) * 1.2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = c("Round 1" = "#2878B5", "Round 2" = "#E07A3F")
  ) +
  labs(
    title = "Multiple antibiotic resistance (MAR) index",
    subtitle = "Mean proportion of tested antibiotics to which isolates were resistant",
    x = "Mean MAR index",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey35"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save MAR plot
save_plot(p8b, "08_ecoli_MAR_by_round", w = 9, h = 4.5)

# Export summary table
write_csv(
  mar_summary %>%
    mutate(
      mean_MAR = round(mean_MAR, 3),
      median_MAR = round(median_MAR, 3)
    ),
  file.path(out_dir, "08_ecoli_MAR_by_round.csv")
)




library(ggplot2)


p8 <- ggplot(
  trend,
  aes(
    x = pct_R,
    y = reorder(antibiotic, pct_R),
    fill = round
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_R)),
    position = position_dodge(width = 0.8),
    hjust = -0.15,
    size = 3.2
  ) +
  scale_fill_manual(
    values = c("Round 1" = "#2C7FB8", "Round 2" = "#D95F0E"),
    name = "Sampling round"
  ) +
  scale_x_continuous(
    labels = label_percent(scale = 1),
    limits = c(0, 110),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "E. coli Antibiotic Resistance by Sampling Round",
    subtitle = "Percentage of tested isolates resistant to each antibiotic",
    x = "Resistant isolates (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )
