library("tidyverse")

library("survival")

library(readxl)

df <- read_excel("D:/Downloads/Homabay_Data.xlsx")


glimpse(df)

Summary(df)

str(df)

Surv(time = df$time, event = df$monthsantibodytest)


table(df$monthsantibodytest, useNA = "ifany")


df$event <- ifelse(df$monthsantibodytest == "positive", 1,
                   ifelse(df$monthsantibodytest == "negative", 0, NA))

df_clean <- df %>%
  filter(!is.na(event))

Surv(time = df$time, event = df$event)



