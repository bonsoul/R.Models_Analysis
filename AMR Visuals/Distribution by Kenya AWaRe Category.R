library(ggplot2)
library(dplyr)
library(scales)

# ---- Data (grand totals by WHO AWaRe group) ----
aware <- tibble(
  group = factor(c("Access", "Watch", "Reserve", "Uncategorized"),
                 levels = c("Uncategorized", "Reserve", "Watch", "Access")),
  n     = c(1503, 1223, 186, 58)
) %>%
  mutate(
    pct         = n / sum(n),
    # label only slices big enough to hold text
    slice_label = ifelse(pct >= 0.05, percent(pct, accuracy = 1), "")
  )

# ---- Colour code (follows the AWaRe traffic-light convention) ----
fill_cols <- c(
  "Access"        = "#2E8B57",   # green  - first-choice, lowest resistance risk
  "Watch"         = "#F2A900",   # amber  - higher resistance potential
  "Reserve"       = "#C0392B",   # red    - last-resort
  "Uncategorized" = "#A6A6A6"    # grey   - not classified
)
text_cols <- c("Access" = "white", "Watch" = "#1A1A1A",
               "Reserve" = "white", "Uncategorized" = "white")

# ---- Plot (donut) ----
p <- ggplot(aware, aes(x = 2, y = n, fill = group)) +
  geom_col(width = 1, colour = "white", linewidth = 2) +
  geom_text(aes(label = slice_label, colour = group),
            position = position_stack(vjust = 0.5),
            fontface = "bold", size = 6, show.legend = FALSE) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +                                           # donut hole
  scale_fill_manual(
    values = fill_cols,
    breaks = c("Access", "Watch", "Reserve", "Uncategorized")
  ) +
  scale_colour_manual(values = text_cols, guide = "none") +
  labs(
    title = "Distribution by Kenya AWaRe Category (2025)",
    fill  = NULL
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title            = element_text(face = "bold", size = 22, colour = "#1A1A1A",
                                         margin = margin(b = 6)),
    plot.title.position   = "plot",
    legend.position       = "right",
    legend.text           = element_text(size = 14, colour = "#333333",
                                         margin = margin(l = 6, b = 6, t = 6)),
    legend.key.size       = unit(0.8, "cm"),
    plot.margin           = margin(25, 25, 15, 25),
    plot.background       = element_rect(fill = "white", colour = NA)
  )

print(p)
# ---- Save ----
#dir.create("outputs", showWarnings = FALSE)
ggsave("01aware_donut.png", p, width = 11, height = 6.5, dpi = 300, bg = "white")
#ggsave("outputs/aware_donut.pdf", p, width = 11, height = 6.5)
