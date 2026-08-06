
##libraries
library(dplyr)
library(tidytuesdayR)
library(tidyverse)

##Importing the datset
nde_experiences <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/main/data/2026/2026-07-21/nde_experiences.csv')


##structure
str(nde_experiences)


## missingness per column
colSums(is.na(nde_experiences)) %>% sort(decreasing = TRUE)


dplyr::glimpse(nde_experiences)


# Greyson score distribution
summary(nde_experiences$greyson_score)
hist(nde_experiences$greyson_score, breaks = 20)

# Classification categories
table(nde_experiences$classification, useNA = "ifany")


# Gender breakdown
table(nde_experiences$gender, useNA = "ifany")


# Country - top 15
nde_experiences %>% count(country, sort = TRUE) %>% print(n = 15)


# Category
table(nde_experiences$category, useNA = "ifany")


# Date ranges
range(nde_experiences$exp_date, na.rm = TRUE)
range(nde_experiences$post_date, na.rm = TRUE)


# Prevalence of each ai_* flag
nde_experiences %>%
  summarise(across(starts_with("ai_"), ~mean(.x, na.rm = TRUE))) %>%
  tidyr::pivot_longer(everything(), names_to = "feature", values_to = "prevalence") %>%
  arrange(desc(prevalence))




library(dplyr)

nde_experiences <- nde_experiences %>%
  mutate(
    classification_clean = case_when(
      classification %in% c("NDE") ~ "NDE",
      classification %in% c("NDE-Like", "NDE-Like;NDE") ~ "NDE-Like",
      classification %in% c("Possible NDE") ~ "Possible NDE",
      classification %in% c("Probable NDE") ~ "Probable NDE",
      classification %in% c("FDE", "SDE", "SMR") ~ "Other",
      TRUE ~ "Other"
    )
  )

table(nde_experiences$classification_clean, useNA = "ifany")


feature_cols <- c("ai_obe", "ai_unity", "ai_hellish",
                  "ai_esp", "ai_past_lives", "ai_world_future", "ai_aliens")

feature_matrix <- nde_experiences %>%
  select(all_of(feature_cols)) %>%
  mutate(across(everything(), as.numeric))

cor_mat <- cor(feature_matrix, use = "pairwise.complete.obs")
round(cor_mat, 2)
