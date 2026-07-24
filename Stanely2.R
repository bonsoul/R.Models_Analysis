library(readxl)
library(tidyverse)
library(gtsummary)

data <- read_excel("C:/Users/LENOVO/Downloads/stanley_data_mapped.xlsx", sheet = "Labeled Data")

att_items <- c("att_a","att_b","att_c","att_d","att_e","att_f")
likert_levels <- c("Strongly Agree","Agree","Neutral","Disagree","Strongly Disagree")

data %>%
  select(all_of(att_items)) %>%
  mutate(across(everything(), ~ factor(.x, levels = likert_levels))) %>%
  tbl_summary()

# 1. Long format + percentages per item
likert_pct <- data %>%
  select(all_of(att_items)) %>%
  pivot_longer(everything(), names_to = "item", values_to = "response") %>%
  filter(!is.na(response), !str_detect(response, "UNRECOGNIZED")) %>%
  mutate(response = factor(response, levels = likert_levels)) %>%
  count(item, response) %>%
  group_by(item) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# 2. Split Neutral in half so it can straddle the zero line
diverging_data <- likert_pct %>%
  mutate(
    pct_signed = case_when(
      response %in% c("Strongly Disagree", "Disagree") ~ -pct,
      response == "Neutral" ~ pct / 2,
      TRUE ~ pct
    )
  )

neutral_split <- diverging_data %>%
  filter(response == "Neutral") %>%
  mutate(side = "left", pct_signed = -pct_signed) %>%
  bind_rows(
    diverging_data %>% filter(response == "Neutral") %>% mutate(side = "right")
  )

diverging_data <- diverging_data %>%
  filter(response != "Neutral") %>%
  mutate(side = NA) %>%
  bind_rows(neutral_split) %>%
  mutate(response = factor(response, levels = likert_levels))

likert_colors <- c(
  "Strongly Disagree" = "#791F1F",
  "Disagree"           = "#E24B4A",
  "Neutral"            = "#B4B2A9",
  "Agree"              = "#639922",
  "Strongly Agree"     = "#27500A"
)

# 4. Plot
ggplot(diverging_data, aes(x = pct_signed, y = item, fill = response)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = likert_colors, name = NULL) +
  scale_x_continuous(labels = function(x) paste0(abs(x), "%")) +
  labs(
    title = "Attitude towards herbal medicine (Q5–Q10)",
    x = "Percent of respondents",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

# --- FIX STARTS HERE: convert text labels back to numeric before scoring ---

score_map <- c("Strongly Agree"=1, "Agree"=2, "Neutral"=3, "Disagree"=4, "Strongly Disagree"=5)

data <- data %>%
  mutate(across(all_of(att_items), ~ recode(.x, !!!score_map), .names = "{.col}_num"))

att_items_num <- paste0(att_items, "_num")

# Valid range check + explicit score computation
att_valid <- data %>%
  select(all_of(att_items_num)) %>%
  mutate(across(everything(), ~ .x %in% 1:5)) %>%
  reduce(`&`)

df <- data %>%
  mutate(
    attitude_score = if_else(
      att_valid,
      rowSums(select(., all_of(att_items_num)), na.rm = FALSE),
      NA_real_
    ),
    attitude_cat = case_when(
      is.na(attitude_score)     ~ NA_character_,
      attitude_score < 18       ~ "Positive",
      attitude_score > 18       ~ "Negative",
      attitude_score == 18      ~ "Borderline",
      TRUE                      ~ NA_character_
    )
  )

# Quick check: how many rows got excluded, and why
sum(is.na(df$attitude_score))
df %>% filter(is.na(attitude_score)) %>% select(all_of(att_items))


library(dplyr)

# Convert Likert responses to numeric
score_map <- c(
  "Strongly Agree" = 1,
  "Agree" = 2,
  "Neutral" = 3,
  "Disagree" = 4,
  "Strongly Disagree" = 5
)

data <- data %>%
  mutate(across(all_of(att_items),
                ~ recode(.x, !!!score_map),
                .names = "{.col}_num"))

att_items_num <- paste0(att_items, "_num")

# Composite score (mean of all attitude items)
df <- data %>%
  mutate(
    attitude_composite = rowMeans(
      select(., all_of(att_items_num)),
      na.rm = FALSE
    )
  )

# View summary
summary(df$attitude_composite)



df <- df %>%
  mutate(
    attitude_cat = case_when(
      attitude_composite < 2.5 ~ "Positive",
      attitude_composite == 2.5 ~ "Borderline",
      attitude_composite > 2.5 ~ "Negative",
      TRUE ~ NA_character_
    )
  )


ggplot(df, aes(x = attitude_composite)) +
  geom_histogram(binwidth = 0.25, fill = "#2C7FB8", color = "white") +
  geom_vline(aes(xintercept = mean(attitude_composite, na.rm = TRUE)),
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = median(attitude_composite, na.rm = TRUE)),
             color = "darkgreen", linetype = "dotted", linewidth = 1) +
  labs(
    title = "Distribution of Composite Attitude Scores",
    subtitle = "Red dashed = Mean; Green dotted = Median",
    x = "Composite Attitude Score",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)
