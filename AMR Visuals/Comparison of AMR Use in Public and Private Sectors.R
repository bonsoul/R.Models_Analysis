library(ggplot2)
library(dplyr)
library(scales)

# ---- Data ----
pps <- tribble(
  ~Year, ~Category, ~Figure_1, ~Figure_2,
  2024, "Private",   125,  257,
  2024, "Public",   1917, 4485,
  2025, "Private",   135,  260,
  2025, "Public",   1850, 4557
)

# ---- Prevalence + 95% CI ----
pps <- pps %>%
  rowwise() %>%
  mutate(
    prev  = Figure_1 / Figure_2 * 100,
    lower = prop.test(Figure_1, Figure_2, correct = FALSE)$conf.int[1] * 100,
    upper = prop.test(Figure_1, Figure_2, correct = FALSE)$conf.int[2] * 100
  ) %>%
  ungroup() %>%
  mutate(
    Year     = factor(Year, levels = c(2024, 2025)),
    # first level is drawn at the bottom of each group, so Private ends up on top
    Category = factor(Category, levels = c("Public", "Private"))
  )

# ---- Plot ----
dodge <- position_dodge(width = 0.7)

p <- ggplot(pps, aes(x = prev, y = Year, fill = Category)) +
  geom_col(position = dodge, width = 0.6) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = dodge, width = 0.2,
                colour = "#444444", linewidth = 0.5) +
  geom_text(aes(x = upper, label = sprintf("%.1f%%", prev)),
            position = dodge, hjust = -0.25,
            size = 4.5, fontface = "bold", colour = "#1A1A1A") +
  scale_fill_manual(values = c("Private" = "#56B4E9", "Public" = "#1F4E79")) +
  scale_y_discrete(limits = rev) +                       # 2024 on top, 2025 below
  scale_x_continuous(limits = c(0, 70),
                     breaks = seq(0, 60, 10),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  guides(fill = guide_legend(reverse = TRUE)) +          # legend reads Private, Public
  labs(
    title    = "Comparison of AMR Use in Public and Private Sectors",
    x = "Prevalence of antimicrobial use", y = NULL, fill = NULL
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
    axis.text.y          = element_text(size = 15, face = "bold", colour = "#333333"),
    axis.text.x          = element_text(colour = "#777777"),
    axis.title.x         = element_text(colour = "#555555", margin = margin(t = 10)),
    legend.position      = "top",
    legend.justification = "left",
    legend.text          = element_text(size = 12),
    plot.caption         = element_text(size = 9, colour = "#888888", hjust = 0,
                                        margin = margin(t = 12)),
    plot.margin          = margin(20, 25, 15, 20),
    plot.background      = element_rect(fill = "white", colour = NA)
  )

print(p)

# ---- Save ----

ggsave("1pps_year_sector_comparison.png", p, width = 9, height = 5.5, dpi = 300, bg = "white")
