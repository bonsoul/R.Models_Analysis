# Total PPS prevalence (all hospitals): previous vs current survey
# install.packages(c("ggplot2", "dplyr"))  # if needed
library(ggplot2)
library(dplyr)

# 1. Data (totals only) --------------------------------------------------
pps <- data.frame(
  survey   = factor(c("Previous PPS", "Current PPS"),
                    levels = c("Previous PPS", "Current PPS")),
  on_abx   = c(2102, 1985),   # patients on antimicrobials
  surveyed = c(4742, 4817)    # patients surveyed
)

# 2. Prevalence + 95% Wilson confidence intervals ------------------------
pps <- pps %>%
  rowwise() %>%
  mutate(
    prev  = on_abx / surveyed * 100,
    lower = prop.test(on_abx, surveyed)$conf.int[1] * 100,
    upper = prop.test(on_abx, surveyed)$conf.int[2] * 100,
    label = sprintf("%.1f%%\n(%s / %s)", prev,
                    format(on_abx, big.mark = ","),
                    format(surveyed, big.mark = ","))
  ) %>%
  ungroup()

# 3. Like-for-like previous total (excluding Mbagathi: 148 / 273) -------
prev_excl <- (2102 - 148) / (4742 - 273) * 100   # 43.7%

# 4. Is the change statistically significant? ----------------------------
print(prop.test(x = pps$on_abx, n = pps$surveyed))   # 2-sample test of proportions
change_pp <- diff(pps$prev)
cat(sprintf("Change: %.1f percentage points (%.1f%% relative)\n",
            change_pp, change_pp / pps$prev[1] * 100))

# 5. Plot -----------------------------------------------------------------
p <- ggplot(pps, aes(x = survey, y = prev, group = 1)) +
  geom_area(fill = "#2a78d6", alpha = 0.10) +
  geom_line(colour = "#2a78d6", linewidth = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.08, colour = "#2a78d6") +
  geom_point(size = 4, colour = "#2a78d6") +
  geom_text(aes(label = label), vjust = -1.6, size = 4, lineheight = 0.95) +
  # like-for-like marker for previous survey
  annotate("point", x = 1, y = prev_excl, shape = 18, size = 4, colour = "grey50") +
  annotate("text", x = 1.05, y = prev_excl - 1.2, hjust = 0, size = 3.3,
           colour = "grey40", label = "43.7% excl. Mbagathi") +
  scale_y_continuous(limits = c(30, 52),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Trend analysis total antimicrobial prevalence",
    subtitle = sprintf("Down %.1f percentage points between surveys",
                       abs(change_pp)),
    x = NULL, y = "Prevalence of antimicrobial use"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold"),
    plot.caption       = element_text(colour = "grey40", hjust = 0)
  )

print(p)

# 6. Save -----------------------------------------------------------------
ggsave("pps_total_trend.png", p, width = 7, height = 5, dpi = 300)