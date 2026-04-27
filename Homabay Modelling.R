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
df1 <- read_excel("D:/Downloads/Homabay_Data.xlsx")





#data manipulation



df <- df1 %>%
  mutate(
    hiv_positive = as.integer(
      birthpcr == "positive" |
        monthspcr == "positive" |
        yearpcr   == "positive" |
        monthsantibodytest == "positive"
    ),
    
    hiv_positive_fac = factor(hiv_positive, levels = c(0, 1),
                              labels = c("HIV Negative", "HIV Positive")),
    
    surv_time = `_t`
  )


df <- df1 %>%
  mutate(
    # ── Outcome ───────────────────────────────────────────────
    hiv_positive = as.integer(
      birthpcr == "positive" |
        monthspcr == "positive" |
        yearpcr   == "positive" |
        monthsantibodytest == "positive"
    ),
    
    hiv_positive_fac = factor(hiv_positive, levels = c(0, 1),
                              labels = c("HIV Negative", "HIV Positive")),
    
    surv_time = `_t`,
    
    # ── Maternal factors ──────────────────────────────────────
    age_group = cut(age, breaks = c(0, 24, 34, Inf),
                    labels = c("≤24 yrs", "25–34 yrs", "≥35 yrs")),
    
    educationlevel = factor(educationlevel,
                            levels = c("never attended school", "primary",
                                       "secondary", "college")),
    
    maritalstatus = factor(maritalstatus,
                           levels = c("Married", "single", "widowed")),
    
    residence = factor(residence, levels = c("rural", "urban")),
    
    employmentstatus = factor(employmentstatus,
                              levels = c("unemployed", "self employed",
                                         "formally employed")),
    
    ancattendance = factor(ancattendance, levels = c("no", "yes")),
    
    haartduringpregnancy = factor(haartduringpregnancy, levels = c("no", "yes")),
    
    adherence = factor(adherence, levels = c("poor", "fair", "good")),
    
    cd4cat = factor(cd4cellcountcellsmm3,
                    levels = c("below 200", "200-500", "above 500")),
    
    whohivdiseasestage = factor(whohivdiseasestage,
                                levels = c("stage I", "stage II",
                                           "stage III", "stage IV")),
    
    hivstatusbeforepregnancy = factor(hivstatusbeforepregnancy,
                                      levels = c("unknown", "positive")),
    
    patnershivstatus = factor(patnershivstatus,
                              levels = c("negative", "unknown", "positive")),
    
    syphillis = factor(syphillis, levels = c("negative", "positive")),
    
    historyofstiduringpregnancy = factor(historyofstiduringpregnancy,
                                         levels = c("no", "yes")),
    
    malaria = factor(malaria, levels = c("no", "yes")),
    
    anaemia = factor(anaemia, levels = c("no", "yes")),
    
    hypertention = factor(hypertention, levels = c("no", "yes")),
    
    treatedduringpregnancy = factor(treatedduringpregnancy, levels = c("no", "yes")),
    
    tmembraner_cat = case_when(
      grepl("spontaneous", tmembraner, ignore.case = TRUE) ~ "Spontaneous ROM",
      grepl("premature",   tmembraner, ignore.case = TRUE) ~ "PROM",
      TRUE ~ "Other"
    ) %>% factor(levels = c("Spontaneous ROM", "PROM", "Other")),
    
    # ── Infant factors ────────────────────────────────────────
    sexofthebaby = factor(sexofthebaby, levels = c("female", "male")),
    
    bwt_cat = cut(bwtkgs, breaks = c(0, 2.5, Inf),
                  labels = c("<2.5 kg (LBW)", "≥2.5 kg"), right = FALSE),
    
    durationofbfmonths_cat = cut(durationofbfmonths,
                                 breaks = c(-1, 6, 12, Inf),
                                 labels = c("≤6 months", "7–12 months", ">12 months")),
    
    ancvisits_cat = cut(numberofancvisitsmade,
                        breaks = c(-1, 0, 3, Inf),
                        labels = c("0 visits", "1–3 visits", "≥4 visits"))
  ) %>%
  
  # ⚠️ FIX: `_st` must use backticks
  filter(`_st` == 1)
    
    
cat("Dataset dimensions:", nrow(df), "rows x", ncol(df), "cols\n")
cat("HIV-positive infants:", sum(df$hiv_positive), "\n")
cat("Censored (HIV-negative):", sum(df$hiv_positive == 0), "\n\n")
    



################################################################################
#  OBJECTIVE 1 – SURVIVAL RATE OF HIV-EXPOSED INFANTS
################################################################################

cat("========== OBJECTIVE 1: SURVIVAL RATE ==========\n")

surv_obj <- Surv(time = df$surv_time, event = df$surv_event)

# ── 1a. Overall Kaplan-Meier ──────────────────────────────────────────────────
km_overall <- survfit(surv_obj ~ 1, data = df)
print(summary(km_overall, times = c(12, 24, 36, 52, 78)))  # weeks

cat("\nOverall survival summary:\n")
cat("Median follow-up time (weeks):", median(df$surv_time), "\n")
cat(sprintf("HIV-free survival rate: %.1f%% (95%% CI: %.1f%% – %.1f%%)\n",
            summary(km_overall)$surv[1] * 100,
            summary(km_overall)$lower[1] * 100,
            summary(km_overall)$upper[1] * 100))

# ── 1b. KM by sub-county ──────────────────────────────────────────────────────
km_subcounty <- survfit(surv_obj ~ subcounty, data = df)

p1a <- ggsurvplot(
  km_overall,
  data          = df,
  risk.table    = TRUE,
  pval          = FALSE,
  conf.int      = TRUE,
  xlab          = "Follow-up Time (Weeks)",
  ylab          = "HIV-Free Survival Probability",
  title         = "Figure 1: Overall HIV-Free Survival of Exposed Infants\n(Homa Bay County)",
  ggtheme       = theme_bw(base_size = 12),
  palette       = "#2C7BB6",
  surv.median.line = "hv",
  risk.table.col   = "strata",
  legend           = "none"
)

p1b <- ggsurvplot(
  km_subcounty,
  data          = df,
  risk.table    = TRUE,
  pval          = TRUE,
  conf.int      = FALSE,
  xlab          = "Follow-up Time (Weeks)",
  ylab          = "HIV-Free Survival Probability",
  title         = "Figure 2: HIV-Free Survival by Sub-County",
  ggtheme       = theme_bw(base_size = 12),
  palette       = "jco",
  legend.title  = "Sub-County",
  legend.labs   = levels(factor(df$subcounty))
)

# Print plots
print(p1a)
print(p1b)

# ── 1c. Cumulative incidence table ────────────────────────────────────────────
ci_table <- summary(km_overall, times = c(4, 8, 12, 24, 52, 78))
ci_df <- data.frame(
  Time_weeks       = ci_table$time,
  N_at_risk        = ci_table$n.risk,
  N_events         = ci_table$n.event,
  Survival         = round(ci_table$surv,   4),
  CI_lower_95      = round(ci_table$lower,  4),
  CI_upper_95      = round(ci_table$upper,  4),
  Cumulative_hazard = round(ci_table$cumhaz, 4)
)
cat("\nTable 1: Kaplan-Meier Survival Estimates at Selected Time Points\n")
print(ci_df)





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



