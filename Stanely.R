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
      response == "Neutral" ~ pct / 2,   # will duplicate below for both sides
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
  bind_rows(neutral_split)

# 3. Set fill order and colors (diverging, matches your color logic)
diverging_data <- diverging_data %>%
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
