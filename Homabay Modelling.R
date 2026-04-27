################################################################################
#  HIV-EXPOSED INFANT SURVIVAL ANALYSIS – HOMA BAY COUNTY
#  Objectives:
#    1. Survival rate of HIV-exposed infants
#    2. Maternal factors associated with HIV antigen positivity
#    3. Infant factors associated with HIV antigen positivity
#    4. HIV-free survival model (Cox PH)
################################################################################



packages <- c("survival", "readxl", "tidyverse",
              "survminer", "finalfit", "gtsummary",
              "flextable", "officer", "car",
              "broom", "scales", "patchwork", "gridExtra")

# Install missing packages
installed <- packages %in% installed.packages()[, "Package"]
if(any(!installed)) install.packages(packages[!installed])

# Load all
lapply(packages, library, character.only = TRUE)




# data
df <- read_excel("D:/Downloads/Homabay_Data.xlsx")


glimpse(df)

summary(df)

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



