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
         category = factor(category, levels = c("Uncategorised", "Reserve", "Watch", "Access")),
         year = factor(year, levels = c(2025, 2024)))  # 2024 on top

ggplot(aware, aes(pct, year, fill = category)) +
  geom_col(width = 0.6, colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5), size = 3.5,
            colour = ifelse(aware$category == "Uncategorised", "black", "white")) +
  scale_fill_manual(values = c("Access" = "#2e7d32", "Watch" = "#f9a825",
                               "Reserve" = "#c62828", "Uncategorised" = "#bdbdbd"),
                    breaks = c("Access", "Watch", "Reserve", "Uncategorised")) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(x = "% of antibiotic prescriptions", y = NULL, fill = NULL,
       title = "AWaRe classification of antibiotics prescribed, 2024 vs 2025",
       caption = "Denominators: 3,198 (2024), 2,970 (2025).") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid = element_blank())

ggsave("12aware_2024_vs_2025.png", width = 9, height = 4, dpi = 300, bg = "white")