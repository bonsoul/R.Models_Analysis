library(tidyverse)

aware <- tribble(
  ~category,        ~year, ~n,   ~total,
  "Access",         2024, 1557, 3198,
  "Watch",          2024, 1385, 3198,
  "Reserve",        2024,  173, 3198,
  "Uncategorised",  2024,   83, 3198,
  "Access",         2025, 1503, 2970,
  "Watch",          2025, 1223, 2970,
  "Reserve",        2025,  186, 2970,
  "Uncategorised",  2025,   58, 2970
) |>
  mutate(pct = n / total * 100,
         category = factor(category, levels = rev(c("Access", "Watch", "Reserve", "Uncategorised"))),
         year = factor(year, levels = c(2025, 2024)))

ggplot(aware, aes(pct, category, fill = year)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_dodge(width = 0.8), hjust = -0.15, size = 3.5) +
  scale_fill_manual(values = c("2024" = "#9ecae1", "2025" = "#08519c"),
                    breaks = c("2024", "2025")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "% of antibiotic prescriptions", y = NULL, fill = NULL,
       title = "AWaRe classification of antibiotics prescribed, 2024 vs 2025",
       caption = "Denominators: 3,198 (2024), 2,970 (2025).") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

ggsave("09aware_2024_vs_2025_hbar.png", width = 8, height = 5, dpi = 300, bg = "white")