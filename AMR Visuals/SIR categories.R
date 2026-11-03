library(tidyverse)

out_dir <- "amr_outputs"
dir.create(out_dir, showWarnings = FALSE)

# Colours for SIR categories
pal_sir <- c(
  S = "#2E8B57",  # Susceptible
  I = "#F2B134",  # Intermediate
  R = "#C0392B"   # Resistant
)

# Default plotting theme
theme_set(
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
)

save_plot <- function(p, name, w = 9, h = 5.5)
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")

# ---- Display the palette ----
pal_sir   # prints the hex codes

palette_df <- tibble(
  code  = factor(names(pal_sir), levels = names(pal_sir)),
  label = c("Susceptible", "Intermediate", "Resistant"),
  hex   = unname(pal_sir)
)

p_pal <- ggplot(palette_df, aes(code, 1, fill = code)) +
  geom_tile(width = 0.95, height = 0.95) +
  geom_text(aes(label = paste0(label, "\n", hex)), colour = "white", fontface = "bold") +
  scale_fill_manual(values = pal_sir, guide = "none") +
  labs(title = "SIR colour palette", x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_pal                          # shows in the Plots pane
# save_plot(p_pal, "00_palette", w = 6, h = 2.5)   # optional: save it