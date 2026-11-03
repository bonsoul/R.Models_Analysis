library(ggplot2)
library(dplyr)
library(scales)

# ---- Data (patients on antimicrobials by age group) ----
age_levels <- c("Adult (≥ 18 years)", "Child (5–17 years)", "Child (< 5 years)",
                "Infant (≤ 12 months)", "Neonate (≤ 28 days)")

age_df <- tibble(
  Year  = rep(c("2024", "2025"), each = 5),
  Age   = rep(age_levels, times = 2),
  n     = c(1313, 190, 166, 162, 271,     # 2024 (N = 2,102)
            1217, 201, 183, 156, 228)     # 2025 (N = 1,985)
) %>%
  group_by(Year) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    Age  = factor(Age, levels = rev(age_levels)),        # Adult on top
    Year = factor(Year, levels = c("2025", "2024"))      # 2024 above 2025 within each group
  )

# ---- Plot ----
dodge <- position_dodge(width = 0.75)

p <- ggplot(age_df, aes(x = pct, y = Age, fill = Year)) +
  geom_col(position = dodge, width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = dodge, hjust = -0.2,
            size = 4.2, fontface = "bold", colour = "#1A1A1A") +
  scale_fill_manual(values = c("2024" = "#9DB8D9", "2025" = "#1F4E79")) +
  scale_x_continuous(
    limits = c(0, 75),
    breaks = seq(0, 75, 15),
    labels = function(x) paste0(x, "%"),
    expand = expansion(add = c(0, 2))
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(
    title    = "Age distribution of patients on antimicrobials",
    subtitle = "2024 (N = 2,102) vs 2025 (N = 1,985)",
    x = "Share of patients on antimicrobials", y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor     = element_blank(),
    panel.grid.major.y   = element_blank(),
    panel.grid.major.x   = element_line(colour = "#E8E8E8"),
    plot.title           = element_text(face = "bold", size = 20, colour = "#1A1A1A",
                                        margin = margin(b = 6)),
    plot.subtitle        = element_text(size = 12, colour = "#555555",
                                        margin = margin(b = 14)),
    axis.text.y          = element_text(size = 13, face = "bold", colour = "#333333"),
    axis.text.x          = element_text(size = 12, colour = "#555555"),
    axis.title.x         = element_text(colour = "#555555", margin = margin(t = 10)),
    legend.position      = "top",
    legend.justification = "left",
    legend.text          = element_text(size = 12),
    plot.caption         = element_text(size = 9, colour = "#888888", hjust = 0,
                                        margin = margin(t = 12)),
    plot.title.position  = "plot",
    plot.caption.position = "plot",
    plot.margin          = margin(20, 25, 15, 20),
    plot.background      = element_rect(fill = "white", colour = NA)
  )

print(p)

# ---- Save ----

ggsave("age_distribution_2024_vs_2025.png", p, width = 9.5, height = 7, dpi = 300, bg = "white")
