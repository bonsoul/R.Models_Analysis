# Sex split of patients on antimicrobials: PPS3 vs PPS4
# install.packages(c("ggplot2", "dplyr"))  # if needed
library(ggplot2)
library(dplyr)

# 1. Data (counts from the totals) ---------------------------------------
sex_df <- data.frame(
  survey = rep(c("PPS3", "PPS4"), each = 2),
  sex    = rep(c("Female", "Male"), times = 2),
  n      = c(1070, 1032,    # PPS3: total 2,102
             1037,  948)    # PPS4: total 1,985
)

# 2. Percentages within each survey --------------------------------------
sex_df <- sex_df %>%
  group_by(survey) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    # PPS3 on top; Male first level so Female sits on the left
    survey = factor(survey, levels = c("PPS4", "PPS3")),
    sex    = factor(sex,    levels = c("Male", "Female"))
  )

print(sex_df)

# 3. Is the change in sex split significant? ------------------------------
female_n <- c(PPS3 = 1070, PPS4 = 1037)
total_n  <- c(PPS3 = 2102, PPS4 = 1985)
print(prop.test(female_n, total_n))

# 4. Plot -------------------------------------------------------------------
p <- ggplot(sex_df, aes(x = pct, y = survey, fill = sex)) +
  geom_col(width = 0.5, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 4.2) +
  scale_fill_manual(values = c(Female = "#6250d6", Male = "#1baf7a"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Sex split of patients on antimicrobials",
    subtitle = "PPS3 (n = 2,102) vs PPS4 (n = 1,985)",
    x = NULL, y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold"),
    legend.position    = "top",
    legend.justification = "left"
  )

print(p)

# 5. Save -------------------------------------------------------------------
ggsave("pps_sex_split.png", p, width = 7, height = 3, dpi = 300)