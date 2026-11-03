
library(ggplot2)
library(dplyr)
library(scales)
library(grid)

# ---- Data ----
kem_df <- tibble(
  kem = factor(c("Yes", "No"), levels = c("Yes", "No")),
  n = c(2914, 56)
) %>%
  mutate(pct = n / sum(n))

total_n <- sum(kem_df$n)

# ---- Slice angles: position No slice on the right ----
rot <- pi / 2 - 2 * pi * (kem_df$pct[1] + kem_df$pct[2] / 2)

kem_df <- kem_df %>%
  mutate(
    end = 2 * pi * cumsum(pct) + rot,
    start = lag(end, default = rot),
    mid = (start + end) / 2,
    explode = c(0, 0.15)
  )

# ---- Construct pie wedges ----
make_wedge <- function(grp, start, end, shift, r = 1, n = 200) {
  t <- seq(start, end, length.out = n)
  mid <- (start + end) / 2
  
  data.frame(
    kem = grp,
    x = c(0, r * sin(t)) + shift * sin(mid),
    y = c(0, r * cos(t)) + shift * cos(mid)
  )
}

slices <- bind_rows(lapply(seq_len(nrow(kem_df)), function(i) {
  make_wedge(
    as.character(kem_df$kem[i]),
    kem_df$start[i],
    kem_df$end[i],
    kem_df$explode[i]
  )
})) %>%
  mutate(kem = factor(kem, levels = levels(kem_df$kem)))

# ---- Label positions ----
yes <- filter(kem_df, kem == "Yes")
no  <- filter(kem_df, kem == "No")

# Main category label outside the pie
yes_label <- data.frame(
  x = -1.18,
  y = 0.10,
  label = paste0(
    "Antibiotic on KEML\n",
    comma(yes$n), " (", percent(yes$pct, accuracy = 0.1), ")"
  )
)

# Callout for the small No slice
dx <- sin(no$mid)
dy <- cos(no$mid)

callout <- data.frame(
  xs = (no$explode + 1.65) * dx,
  ys = (no$explode + 1.65) * dy,
  xe = (no$explode + 1.02) * dx,
  ye = (no$explode + 1.02) * dy
)

blue   <- "#1F4E79"
orange <- "#E07B39"

# ---- Plot ----
p <- ggplot() +
  geom_polygon(
    data = slices,
    aes(x = x, y = y, fill = kem, group = kem),
    colour = "white",
    linewidth = 1.3
  ) +
  
  # Yes category label outside the pie
  annotate(
    "text",
    x = yes_label$x,
    y = yes_label$y,
    label = yes_label$label,
    hjust = 1,
    size = 4.5,
    lineheight = 1.2,
    fontface = "bold",
    colour = blue
  ) +
  
  # No slice callout
  geom_segment(
    data = callout,
    aes(x = xs, y = ys, xend = xe, yend = ye),
    colour = orange,
    linewidth = 0.9,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  
  annotate(
    "text",
    x = callout$xs + 0.08,
    y = callout$ys + 0.10,
    label = "Antibiotic not on KEML",
    hjust = 0,
    size = 4.2,
    fontface = "bold",
    colour = orange
  ) +
  
  annotate(
    "text",
    x = callout$xs + 0.08,
    y = callout$ys - 0.10,
    label = paste0(
      comma(no$n), " (", percent(no$pct, accuracy = 0.1), ")"
    ),
    hjust = 0,
    size = 4,
    colour = orange
  ) +
  
  scale_fill_manual(
    values = c("Yes" = blue, "No" = orange),
    labels = c(
      "Yes" = "On KEML",
      "No" = "Not on KEML"
    )
  ) +
  
  coord_fixed(
    xlim = c(-2.0, 3.0),
    ylim = c(-1.35, 1.35),
    clip = "off"
  ) +
  
  labs(
    title = "Antibiotic Prescriptions by KEML Listing (N = 2,970)",
    fill = NULL
  ) +
  
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold", size = 19, colour = "#1A1A1A",
      margin = margin(b = 7)
    ),
    plot.subtitle = element_text(
      size = 11.5, colour = "#555555",
      margin = margin(b = 18)
    ),
    plot.caption = element_text(
      size = 9, colour = "#777777", hjust = 0,
      margin = margin(t = 12)
    ),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    legend.position = "bottom",
    legend.text = element_text(
      size = 11, colour = "#333333"
    ),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(20, 35, 15, 35),
    plot.background = element_rect(
      fill = "white", colour = NA
    )
  )

print(p)



ggsave(
  filename = "Antibiotics_KEML.png",
  plot = p,
  width = 12,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

