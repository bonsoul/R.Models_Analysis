library(ggplot2)
library(dplyr)
library(scales)

# ---- Data ----
total_df <- tibble(
  ordered = factor(c("No", "Yes"), levels = c("No", "Yes")),
  n       = c(1597, 388)
) %>%
  mutate(
    pct   = n / sum(n),
    label = paste0(percent(pct, accuracy = 1), "\n", comma(n))
  )

total_n <- sum(total_df$n)  # 1,985

# ---- Plot ----
p <- ggplot(total_df, aes(x = "", y = n, fill = ordered)) +
  geom_col(width = 1, colour = "white", linewidth = 1.5) +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 7, lineheight = 0.9) +
  scale_fill_manual(
    values = c("Yes" = "#1F4E79", "No" = "#E07B39"),
    labels = c("No"  = "No culture / AST ordered",
               "Yes" = "Culture / AST ordered")
  ) +
  labs(
    title    = "Culture test and AST ordered (N = 1985)",
    fill     = NULL
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 20, colour = "#1A1A1A",
                                    margin = margin(b = 6)),
    plot.subtitle    = element_text(size = 12, colour = "#555555",
                                    margin = margin(b = 10)),
    plot.caption     = element_text(size = 9, colour = "#888888", hjust = 0,
                                    margin = margin(t = 10)),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12, colour = "#333333"),
    legend.key.size  = unit(0.6, "cm"),
    plot.margin      = margin(20, 20, 15, 20),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

p

# ---- Save ----
#dir.create("outputs", showWarnings = FALSE)

# High-res PNG (slides, email, reports)
ggsave("1culture_ast_total_pie.png", plot = p,
       width = 8, height = 7, dpi = 300, bg = "white")

# Vector formats (crisp at any size, good for print / PowerPoint)
#ggsave("outputs/culture_ast_total_pie.pdf", plot = p,width = 8, height = 7)
#ggsave("outputs/culture_ast_total_pie.svg", plot = p,width = 8, height = 7)   # needs the svglite package