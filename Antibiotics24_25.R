library(tidyverse)

df <- tribble(
  ~drug,                                   ~year, ~n,  ~total,
  "Ceftriaxone",                           2024, 762, 3198,
  "Metronidazole",                         2024, 515, 3198,
  "Amoxicillin-clavulanate",               2024, 256, 3198,
  "Gentamicin",                            2024, 245, 3198,
  "Benzylpenicillin",                      2024, 224, 3198,
  "Meropenem",                             2024, 139, 3198,
  "Flucloxacillin",                        2024, 120, 3198,
  "Clindamycin",                           2024, 113, 3198,
  "Ceftazidime",                           2024, 103, 3198,
  "Cefuroxime",                            2024,  98, 3198,
  "Ceftriaxone",                           2025, 602, 2970,
  "Metronidazole",                         2025, 459, 2970,
  "Amoxicillin-clavulanate",               2025, 388, 2970,
  "Benzylpenicillin",                      2025, 181, 2970,
  "Gentamicin",                            2025, 181, 2970,
  "Meropenem",                             2025, 123, 2970,
  "Ceftazidime",                           2025, 108, 2970,
  "Piperacillin-enzyme inhibitor",         2025, 102, 2970,
  "Flucloxacillin",                        2025, 100, 2970,
  "Clindamycin",                           2025,  91, 2970
) |>
  mutate(pct = n / total * 100,
         year = factor(year)) |>
  complete(drug, year)   # fills missing combos with NA

ord <- df |> group_by(drug) |> summarise(m = mean(pct, na.rm = TRUE)) |>
  arrange(m) |> pull(drug)

library(tidyverse)

ggplot(df, aes(pct, factor(drug, levels = ord), fill = year)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    na.rm = TRUE
  ) +
  geom_text(
    aes(label = ifelse(is.na(pct), "", sprintf("%.1f%%", pct))),
    position = position_dodge(width = 0.8),
    hjust = -0.1,
    size = 3
  ) +
  scale_fill_manual(
    values = c("2024" = "#9ecae1", "2025" = "#08519c")
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "% of antibiotic prescriptions",
    y = NULL,
    fill = NULL,
    title = "Top 10 antibiotics prescribed, 2024 vs 2025"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )


ggsave("01antibiotics_2024_vs_2025.png", width = 10, height = 6, dpi = 300, bg = "white")


