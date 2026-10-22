library(ggplot2)
library(dplyr)
library(scales)

# 1. Data (totals only) ---------------------------------------------------
pps <- tibble(
  survey   = factor(c("Previous PPS\n(2024)", "Current PPS\n(2025)"),
                    levels = c("Previous PPS\n(2024)", "Current PPS\n(2025)")),
  on_abx   = c(2102, 1985),   # patients on antimicrobials
  surveyed = c(4742, 4817)    # patients surveyed
)

# 2. Prevalence + 95% Wilson confidence intervals -------------------------
pps <- pps %>%
  rowwise() %>%
  mutate(
    prev  = on_abx / surveyed * 100,
    lower = prop.test(on_abx, surveyed, correct = FALSE)$conf.int[1] * 100,
    upper = prop.test(on_abx, surveyed, correct = FALSE)$conf.int[2] * 100,
    label = sprintf("%.1f%%\n(%s / %s)", prev,
                    comma(on_abx), comma(surveyed))
  ) %>%
  ungroup()

# 3. Is the change statistically significant? ------------------------------
test      <- prop.test(x = pps$on_abx, n = pps$surveyed)
change_pp <- diff(pps$prev)
p_txt     <- ifelse(test$p.value < 0.001, "p < 0.001",
                    sprintf("p = %.3f", test$p.value))
print(test)
cat(sprintf("Change: %.1f percentage points (%.1f%% relative)\n",
            change_pp, change_pp / pps$prev[1] * 100))

# 4. Plot --------------------------------------------------------------------
blue <- "#1F4E79"

p <- ggplot(pps, aes(x = survey, y = prev, group = 1)) +
  geom_area(fill = blue, alpha = 0.08) +
  geom_line(colour = blue, linewidth = 1.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.07, colour = blue, linewidth = 0.8) +
  geom_point(size = 5, colour = blue, fill = "white", shape = 21, stroke = 1.8) +
  geom_text(aes(y = upper, label = label),
            vjust = -0.6, size = 4.5, fontface = "bold",
            colour = "#1A1A1A", lineheight = 0.95) +
  scale_y_continuous(limits = c(30, 52),
                     breaks = seq(30, 50, 5),
                     labels = function(x) paste0(x, "%")) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  labs(
    title    = "Prevalence of Antimicrobial use ",
    x = NULL, y = "Prevalence of antimicrobial use",
    caption  = "Source: your data source here"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#E8E8E8"),
    plot.title         = element_text(face = "bold", size = 20, colour = "#1A1A1A",
                                      margin = margin(b = 6)),
    plot.subtitle      = element_text(size = 12, colour = "#555555",
                                      margin = margin(b = 14)),
    axis.text.x        = element_text(size = 13, face = "bold", colour = "#333333"),
    axis.text.y        = element_text(colour = "#777777"),
    axis.title.y       = element_text(colour = "#555555", margin = margin(r = 10)),
    plot.caption       = element_text(size = 9, colour = "#888888", hjust = 0,
                                      margin = margin(t = 12)),
    plot.margin        = margin(20, 25, 15, 20),
    plot.background    = element_rect(fill = "white", colour = NA)
  )

print(p)

# 5. Save --------------------------------------------------------------------


ggsave("pps_total_trend1.png", p, width = 8, height = 6, dpi = 300, bg = "white")
#ggsave("outputs/pps_total_trend.pdf", p, width = 8, height = 6)