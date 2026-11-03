library(ggplot2)
library(dplyr)
library(scales)

# ---- Data ----
kem_df <- tibble(
  kem = factor(c("Yes", "No"), levels = c("Yes", "No")),
  n   = c(2914, 56)
) %>%
  mutate(pct = n / sum(n))

total_n <- sum(kem_df$n)                                   # 2,970

# ---- Angles: rotate so the No slice points to the right ----
rot <- pi / 2 - 2 * pi * (kem_df$pct[1] + kem_df$pct[2] / 2)

kem_df <- kem_df %>%
  mutate(
    end     = 2 * pi * cumsum(pct) + rot,
    start   = dplyr::lag(end, default = rot),
    mid     = (start + end) / 2,
    explode = c(0, 0.15)                                   # how far No is pulled out
  )

# ---- Build each wedge as a polygon (no extra packages) ----
make_wedge <- function(grp, start, end, shift, r = 1, n = 200) {
  t <- seq(start, end, length.out = n)
  mid <- (start + end) / 2
  data.frame(
    kem = grp,
    x   = c(0, r * sin(t)) + shift * sin(mid),
    y   = c(0, r * cos(t)) + shift * cos(mid)
  )
}

slices <- bind_rows(lapply(seq_len(nrow(kem_df)), function(i) {
  make_wedge(as.character(kem_df$kem[i]), kem_df$start[i],
             kem_df$end[i], kem_df$explode[i])
})) %>%
  mutate(kem = factor(kem, levels = levels(kem_df$kem)))

# ---- Label positions ----
yes <- kem_df[kem_df$kem == "Yes", ]
no  <- kem_df[kem_df$kem == "No", ]
dx  <- sin(no$mid); dy <- cos(no$mid)                      # direction of the No slice

callout <- data.frame(
  xs = (no$explode + 1.60) * dx, ys = (no$explode + 1.60) * dy,   # arrow start (label side)
  xe = (no$explode + 1.06) * dx, ye = (no$explode + 1.06) * dy    # arrow tip (slice edge)
)

blue   <- "#1F4E79"
orange <- "#E07B39"

# ---- Plot ----
p <- ggplot() +
  geom_polygon(data = slices, aes(x = x, y = y, fill = kem, group = kem),
               colour = "white", linewidth = 1.5) +
  # Yes label in the centre
  annotate("text", x = 0, y = 0.10, label = percent(yes$pct, accuracy = 0.1),
           size = 16, fontface = "bold", colour = "white") +
  annotate("text", x = 0, y = -0.20, label = paste0("Yes  ", comma(yes$n)),
           size = 5.5, colour = "white") +
  # No label outside, with arrow
  geom_segment(data = callout, aes(x = xs, y = ys, xend = xe, yend = ye),
               colour = orange, linewidth = 0.8,
               arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
  annotate("text", x = callout$xs + 0.08, y = callout$ys + 0.10, label = "No",
           hjust = 0, size = 6, fontface = "bold", colour = orange) +
  annotate("text", x = callout$xs + 0.08, y = callout$ys - 0.10,
           label = paste0(comma(no$n), " (", percent(no$pct, accuracy = 0.1), ")"),
           hjust = 0, size = 5, colour = orange) +
  scale_fill_manual(
    values = c("Yes" = blue, "No" = orange),
    labels = c("Yes" = "Antibiotic on KEML", "No" = "Antibiotic not on KEML")
  ) +
  coord_fixed(xlim = c(-1.15, 2.7), ylim = c(-1.2, 1.2), clip = "off") +
  labs(
    title    = "98% of antibiotics are on the Kenya Essential Medicines List",
    subtitle = paste0("Antibiotics by KEML listing (N = ", comma(total_n), ")"),
    fill     = NULL,
    caption  = "Source: your data source here"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title            = element_text(face = "bold", size = 20, colour = "#1A1A1A",
                                         margin = margin(b = 6)),
    plot.subtitle         = element_text(size = 12, colour = "#555555",
                                         margin = margin(b = 4)),
    plot.caption          = element_text(size = 9, colour = "#888888", hjust = 0,
                                         margin = margin(t = 10)),
    plot.title.position   = "plot",
    plot.caption.position = "plot",
    legend.position       = "bottom",
    legend.text           = element_text(size = 13, colour = "#333333",
                                         margin = margin(r = 14)),
    legend.key.size       = unit(0.6, "cm"),
    plot.margin           = margin(25, 25, 15, 25),
    plot.background       = element_rect(fill = "white", colour = NA)
  )

print(p)

# ---- Save ----
dir.create("outputs", showWarnings = FALSE)
ggsave("outputs/keml_pie.png", p, width = 10, height = 6.5, dpi = 300, bg = "white")
ggsave("outputs/keml_pie.pdf", p, width = 10, height = 6.5)