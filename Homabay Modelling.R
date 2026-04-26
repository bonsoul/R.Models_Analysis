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


survival_model <- survfit(Surv(time, event) ~ 1, data = df_clean)


summary(survival_model)


plot(survival_model,
     xlab = "Time",
     ylab = "Survival Probability",
     main = "Kaplan-Meier Curve",
     col = "blue",
     lwd = 2)


survival_model <- survfit(Surv(time, event) ~ patnershivstatus, data = df_clean)

plot(survival_model,
     col = c("red", "blue"),
     lwd = 2,
     xlab = "Time",
     ylab = "Survival Probability")

legend("bottomleft",
       legend = levels(df_clean$patnershivstatus),
       col = c("red", "blue"),
       lwd = 2)
