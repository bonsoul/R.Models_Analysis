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
              "broom", "scales", "patchwork", "gridExtra","sf", "ggplot2", "dplyr", "readr", "readxl", "RColorBrewer",
              "viridis", "patchwork", "ggspatial", "leaflet", "leaflet.extras",
              "htmlwidgets", "scales", "ggrepel", "tidyr", "KernSmooth",
              "rnaturalearth", "rnaturalearthdata")




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
    




df$event <- ifelse(df$monthsantibodytest == "positive", 1,
                   ifelse(df$monthsantibodytest == "negative", 0, NA))

df <- df %>%
  filter(!is.na(event), !is.na(time))



################################################################################
#  OBJECTIVE 1 – SURVIVAL RATE OF HIV-EXPOSED INFANTS
################################################################################

cat("========== OBJECTIVE 1: SURVIVAL RATE ==========\n")

surv_obj <- Surv(time = df$surv_time, event = df$event)

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



################################################################################
#  OBJECTIVE 2 – MATERNAL FACTORS ASSOCIATED WITH HIV ANTIGEN POSITIVITY
################################################################################

cat("\n========== OBJECTIVE 2: MATERNAL FACTORS ==========\n")

# Variables list
maternal_vars <- c(
  "age_group", "educationlevel", "maritalstatus", "residence",
  "employmentstatus", "ancattendance", "ancvisits_cat",
  "haartduringpregnancy", "adherence", "cd4cat", "whohivdiseasestage",
  "hivstatusbeforepregnancy", "patnershivstatus",
  "syphillis", "historyofstiduringpregnancy",
  "malaria", "anaemia", "hypertention", "treatedduringpregnancy",
  "tmembraner_cat"
)

# ── 2a. Descriptive table (Table 1) ──────────────────────────────────────────
tbl_maternal <- df %>%
  select(all_of(c(maternal_vars, "hiv_positive_fac"))) %>%
  tbl_summary(
    by        = hiv_positive_fac,
    statistic = list(all_categorical() ~ "{n} ({p}%)",
                     all_continuous()  ~ "{mean} ({sd})"),
    missing   = "no",
    label     = list(
      age_group                   ~ "Age Group",
      educationlevel              ~ "Education Level",
      maritalstatus               ~ "Marital Status",
      residence                   ~ "Residence",
      employmentstatus            ~ "Employment Status",
      ancattendance               ~ "ANC Attendance",
      ancvisits_cat               ~ "Number of ANC Visits",
      haartduringpregnancy        ~ "HAART During Pregnancy",
      adherence                   ~ "Adherence to HAART",
      cd4cat                      ~ "CD4 Count (cells/mm³)",
      whohivdiseasestage          ~ "WHO HIV Disease Stage",
      hivstatusbeforepregnancy    ~ "HIV Status Before Pregnancy",
      patnershivstatus            ~ "Partner HIV Status",
      syphillis                   ~ "Syphilis",
      historyofstiduringpregnancy ~ "History of STI During Pregnancy",
      malaria                     ~ "Malaria",
      anaemia                     ~ "Anaemia",
      hypertention                ~ "Hypertension",
      treatedduringpregnancy      ~ "Treated for STI During Pregnancy",
      tmembraner_cat              ~ "Rupture of Membranes"
    )
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test",
                    all_continuous()  ~ "t.test")) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 2: Maternal Factors by Infant HIV Status**")

print(tbl_maternal)

# ── 2b. Univariable logistic regression (maternal) ───────────────────────────
univ_maternal <- df %>%
  select(all_of(c(maternal_vars, "hiv_positive"))) %>%
  tbl_uvregression(
    method    = glm,
    y         = hiv_positive,
    method.args = list(family = binomial),
    exponentiate = TRUE,
    label     = list(
      age_group                   ~ "Age Group",
      educationlevel              ~ "Education Level",
      maritalstatus               ~ "Marital Status",
      residence                   ~ "Residence",
      employmentstatus            ~ "Employment Status",
      ancattendance               ~ "ANC Attendance",
      ancvisits_cat               ~ "Number of ANC Visits",
      haartduringpregnancy        ~ "HAART During Pregnancy",
      adherence                   ~ "Adherence to HAART",
      cd4cat                      ~ "CD4 Count (cells/mm³)",
      whohivdiseasestage          ~ "WHO HIV Disease Stage",
      hivstatusbeforepregnancy    ~ "HIV Status Before Pregnancy",
      patnershivstatus            ~ "Partner HIV Status",
      syphillis                   ~ "Syphilis",
      historyofstiduringpregnancy ~ "History of STI During Pregnancy",
      malaria                     ~ "Malaria",
      anaemia                     ~ "Anaemia",
      hypertention                ~ "Hypertension",
      treatedduringpregnancy      ~ "Treated for STI During Pregnancy",
      tmembraner_cat              ~ "Rupture of Membranes"
    )
  ) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Table 3: Univariable Logistic Regression – Maternal Factors**")

print(univ_maternal)

# ── 2c. Multivariable logistic regression (maternal) ─────────────────────────
# Select variables significant at p < 0.25 in univariable (common threshold)
# Using the key clinical & significant variables
mv_maternal_formula <- hiv_positive ~ adherence + cd4cat + whohivdiseasestage +
  haartduringpregnancy + ancattendance + historyofstiduringpregnancy +
  syphillis + patnershivstatus + tmembraner_cat

mv_maternal_model <- glm(mv_maternal_formula, data = df, family = binomial)
cat("\nMultivariable Logistic Regression – Maternal Factors\n")
print(summary(mv_maternal_model))

tbl_mv_maternal <- tbl_regression(
  mv_maternal_model,
  exponentiate = TRUE,
  label = list(
    adherence                   ~ "Adherence to HAART",
    cd4cat                      ~ "CD4 Count (cells/mm³)",
    whohivdiseasestage          ~ "WHO HIV Disease Stage",
    haartduringpregnancy        ~ "HAART During Pregnancy",
    ancattendance               ~ "ANC Attendance",
    historyofstiduringpregnancy ~ "History of STI During Pregnancy",
    syphillis                   ~ "Syphilis",
    patnershivstatus            ~ "Partner HIV Status",
    tmembraner_cat              ~ "Rupture of Membranes"
  )
) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Table 4: Multivariable Logistic Regression – Maternal Factors**")

print(tbl_mv_maternal)




plot(survival_model,
     xlab = "Time",
     ylab = "Survival Probability",
     main = "Kaplan-Meier Curve",
     col = "blue",
     lwd = 2)



################################################################################
#  OBJECTIVE 3 – INFANT FACTORS ASSOCIATED WITH HIV ANTIGEN POSITIVITY
################################################################################

cat("\n========== OBJECTIVE 3: INFANT FACTORS ==========\n")

infant_vars <- c(
  "sexofthebaby", "bwt_cat", "durationofbfmonths_cat",
  "babywt6wks"  # weight at 6 weeks as proxy of growth/nutrition
)

# ── 3a. Descriptive table ────────────────────────────────────────────────────
tbl_infant_desc <- df %>%
  select(all_of(c(infant_vars, "hiv_positive_fac"))) %>%
  tbl_summary(
    by        = hiv_positive_fac,
    statistic = list(all_categorical() ~ "{n} ({p}%)",
                     all_continuous()  ~ "{mean} ({sd})"),
    missing   = "no",
    label     = list(
      sexofthebaby            ~ "Sex of Infant",
      bwt_cat                 ~ "Birth Weight Category",
      durationofbfmonths_cat  ~ "Duration of Breastfeeding",
      babywt6wks              ~ "Infant Weight at 6 Weeks (kg)"
    )
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test",
                    all_continuous()  ~ "t.test")) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 5: Infant Factors by HIV Status**")

print(tbl_infant_desc)

# ── 3b. Univariable logistic regression (infant) ─────────────────────────────
univ_infant <- df %>%
  select(all_of(c(infant_vars, "hiv_positive"))) %>%
  tbl_uvregression(
    method       = glm,
    y            = hiv_positive,
    method.args  = list(family = binomial),
    exponentiate = TRUE,
    label        = list(
      sexofthebaby            ~ "Sex of Infant",
      bwt_cat                 ~ "Birth Weight Category",
      durationofbfmonths_cat  ~ "Duration of Breastfeeding",
      babywt6wks              ~ "Infant Weight at 6 Weeks (kg)"
    )
  ) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Table 6: Univariable Logistic Regression – Infant Factors**")

print(univ_infant)

# ── 3c. Multivariable logistic regression (infant) ────────────────────────────
mv_infant_formula <- hiv_positive ~ sexofthebaby + bwt_cat +
  durationofbfmonths_cat + babywt6wks

mv_infant_model <- glm(mv_infant_formula, data = df, family = binomial)
cat("\nMultivariable Logistic Regression – Infant Factors\n")
print(summary(mv_infant_model))

tbl_mv_infant <- tbl_regression(
  mv_infant_model,
  exponentiate = TRUE,
  label = list(
    sexofthebaby            ~ "Sex of Infant",
    bwt_cat                 ~ "Birth Weight Category",
    durationofbfmonths_cat  ~ "Duration of Breastfeeding",
    babywt6wks              ~ "Infant Weight at 6 Weeks (kg)"
  )
) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Table 7: Multivariable Logistic Regression – Infant Factors**")

print(tbl_mv_infant)

# ── 3d. Forest plot – Infant OR ──────────────────────────────────────────────
mv_infant_tidy <- broom::tidy(mv_infant_model, exponentiate = TRUE,
                              conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(significant = ifelse(p.value < 0.05, "Significant (p<0.05)", "Not significant"))

p3 <- ggplot(mv_infant_tidy, aes(x = estimate, y = reorder(term, estimate),
                                 color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10(labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = c("Significant (p<0.05)" = "#D73027",
                                "Not significant" = "#4393C3")) +
  labs(title = "Figure 4: Forest Plot – Infant Factors (Multivariable OR)",
       x = "Odds Ratio (log scale)",
       y = "",
       color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

print(p3)

# ── 3e. Bar chart: HIV positivity by breastfeeding duration ──────────────────
bf_summary <- df %>%
  filter(!is.na(durationofbfmonths_cat)) %>%
  group_by(durationofbfmonths_cat, hiv_positive_fac) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(durationofbfmonths_cat) %>%
  mutate(pct = n / sum(n) * 100)

p3b <- ggplot(bf_summary, aes(x = durationofbfmonths_cat, y = pct,
                              fill = hiv_positive_fac)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("HIV Negative" = "#4393C3",
                               "HIV Positive" = "#D73027")) +
  labs(title = "Figure 5: Infant HIV Status by Breastfeeding Duration",
       x = "Duration of Breastfeeding",
       y = "Percentage (%)",
       fill = "HIV Status") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p3b)



################################################################################
#  OBJECTIVE 4 – HIV-FREE SURVIVAL MODEL (COX PROPORTIONAL HAZARDS)
################################################################################

cat("\n========== OBJECTIVE 4: HIV-FREE SURVIVAL MODEL ==========\n")

# ── 4a. Test PH assumption before fitting ─────────────────────────────────────
cat("\n--- Schoenfeld Residual Test (PH Assumption) ---\n")
# Fit simple Cox first for PH test
cox_test_vars <- c("adherence", "cd4cat", "whohivdiseasestage",
                   "haartduringpregnancy", "durationofbfmonths_cat",
                   "bwt_cat", "sexofthebaby", "syphillis", "tmembraner_cat")

cox_ph_test_formula <- as.formula(
  paste("surv_obj ~", paste(cox_test_vars, collapse = " + "))
)

cox_ph_test_model <- coxph(cox_ph_test_formula, data = df)
ph_test <- cox.zph(cox_ph_test_model)
print(ph_test)
cat("\nIf p > 0.05 for each variable, PH assumption holds.\n")

# ── 4b. Univariable Cox regression ────────────────────────────────────────────
all_cox_vars <- c(maternal_vars, infant_vars)

tbl_cox_univ <- df %>%
  select(all_of(c(all_cox_vars, "surv_time", "event"))) %>%
  tbl_uvregression(
    method       = coxph,
    y            = Surv(surv_time, event),
    exponentiate = TRUE,
    label        = list(
      age_group                   ~ "Age Group",
      educationlevel              ~ "Education Level",
      maritalstatus               ~ "Marital Status",
      residence                   ~ "Residence",
      employmentstatus            ~ "Employment Status",
      ancattendance               ~ "ANC Attendance",
      ancvisits_cat               ~ "Number of ANC Visits",
      haartduringpregnancy        ~ "HAART During Pregnancy",
      adherence                   ~ "Adherence to HAART",
      cd4cat                      ~ "CD4 Count (cells/mm³)",
      whohivdiseasestage          ~ "WHO HIV Disease Stage",
      hivstatusbeforepregnancy    ~ "HIV Status Before Pregnancy",
      patnershivstatus            ~ "Partner HIV Status",
      syphillis                   ~ "Syphilis",
      historyofstiduringpregnancy ~ "History of STI During Pregnancy",
      malaria                     ~ "Malaria",
      anaemia                     ~ "Anaemia",
      hypertention                ~ "Hypertension",
      treatedduringpregnancy      ~ "Treated for STI During Pregnancy",
      tmembraner_cat              ~ "Rupture of Membranes",
      sexofthebaby                ~ "Sex of Infant",
      bwt_cat                     ~ "Birth Weight Category",
      durationofbfmonths_cat      ~ "Duration of Breastfeeding",
      babywt6wks                  ~ "Infant Weight at 6 Weeks (kg)"
    )
  ) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Table 8: Univariable Cox Regression – All Factors**")

print(tbl_cox_univ)

# ── 4c. Final multivariable Cox PH model ──────────────────────────────────────
# Include variables with clinical plausibility & p < 0.25 from univariable Cox
cox_final_formula <- Surv(surv_time, event) ~
  adherence + cd4cat + whohivdiseasestage + haartduringpregnancy +
  syphillis + historyofstiduringpregnancy + patnershivstatus +
  durationofbfmonths_cat + bwt_cat + sexofthebaby + tmembraner_cat

cox_final_model <- coxph(cox_final_formula, data = df, ties = "efron")
cat("\n--- Final Cox PH Model Summary ---\n")
print(summary(cox_final_model))

# ── 4d. gtsummary table for final Cox model ───────────────────────────────────
tbl_cox_final <- tbl_regression(
  cox_final_model,
  exponentiate = TRUE,
  label = list(
    adherence                   ~ "Adherence to HAART",
    cd4cat                      ~ "CD4 Count (cells/mm³)",
    whohivdiseasestage          ~ "WHO HIV Disease Stage",
    haartduringpregnancy        ~ "HAART During Pregnancy",
    syphillis                   ~ "Syphilis",
    historyofstiduringpregnancy ~ "History of STI During Pregnancy",
    patnershivstatus            ~ "Partner HIV Status",
    durationofbfmonths_cat      ~ "Duration of Breastfeeding",
    bwt_cat                     ~ "Birth Weight Category",
    sexofthebaby                ~ "Sex of Infant",
    tmembraner_cat              ~ "Rupture of Membranes"
  )
) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Table 9: Final Multivariable Cox PH Model – HIV-Free Survival**")

print(tbl_cox_final)

# ── 4e. Model diagnostics ─────────────────────────────────────────────────────
cat("\n--- Cox PH Diagnostics ---\n")

# Schoenfeld residuals for final model
ph_final <- cox.zph(cox_final_model)
cat("Global PH test p-value:", ph_final$table["GLOBAL", "p"], "\n")
print(ph_final)

# Plot Schoenfeld residuals
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
plot(ph_final)
par(mfrow = c(1, 1))

# ── 4f. KM plots stratified by key predictors ─────────────────────────────────
km_adherence <- survfit(surv_obj ~ adherence, data = df)
p4a <- ggsurvplot(
  km_adherence, data = df,
  pval          = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  xlab          = "Follow-up Time (Weeks)",
  ylab          = "HIV-Free Survival Probability",
  title         = "Figure 6: HIV-Free Survival by HAART Adherence",
  legend.title  = "Adherence",
  palette       = c("#D73027", "#FEE090", "#4393C3"),
  ggtheme       = theme_bw(base_size = 12)
)
print(p4a)

km_bf <- survfit(surv_obj ~ durationofbfmonths_cat, data = df)
p4b <- ggsurvplot(
  km_bf, data = df,
  pval          = TRUE,
  conf.int      = FALSE,
  risk.table    = TRUE,
  xlab          = "Follow-up Time (Weeks)",
  ylab          = "HIV-Free Survival Probability",
  title         = "Figure 7: HIV-Free Survival by Duration of Breastfeeding",
  legend.title  = "Breastfeeding Duration",
  palette       = "jco",
  ggtheme       = theme_bw(base_size = 12)
)
print(p4b)


km_cd4 <- survfit(surv_obj ~ cd4cat, data = df)
p4c <- ggsurvplot(
  km_cd4, data = df,
  pval          = TRUE,
  conf.int      = FALSE,
  risk.table    = TRUE,
  xlab          = "Follow-up Time (Weeks)",
  ylab          = "HIV-Free Survival Probability",
  title         = "Figure 8: HIV-Free Survival by Maternal CD4 Count",
  legend.title  = "CD4 Count",
  palette       = "npg",
  ggtheme       = theme_bw(base_size = 12)
)
print(p4c)

# ── 4g. Forest plot for final Cox model ──────────────────────────────────────
cox_tidy <- broom::tidy(cox_final_model, exponentiate = TRUE,
                        conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(significant = ifelse(p.value < 0.05, "Significant (p<0.05)", "Not significant"))

p4d <- ggplot(cox_tidy, aes(x = estimate, y = reorder(term, estimate),
                            color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10(labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = c("Significant (p<0.05)" = "#D73027",
                                "Not significant" = "#4393C3")) +
  labs(title = "Figure 9: Forest Plot – Final Cox PH Model\n(HIV-Free Survival)",
       x = "Hazard Ratio (log scale, 95% CI)",
       y = "",
       color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

print(p4d)






################################################################################
# COMBINED MERGED TABLE (Univariable + Multivariable)
################################################################################

cat("\n========== MERGED TABLES ==========\n")

# Merged maternal logistic regression table
tbl_maternal_merged <- tbl_merge(
  tbls        = list(univ_maternal, tbl_mv_maternal),
  tab_spanner = c("**Univariable**", "**Multivariable**")
) %>%
  modify_caption("**Table 10: Logistic Regression – Maternal Factors (Combined)**")

print(tbl_maternal_merged)

# Merged Cox table (univariable + multivariable subset)
tbl_cox_mv_sub <- df %>%
  select(surv_time, event,
         adherence, cd4cat, whohivdiseasestage, haartduringpregnancy,
         syphillis, historyofstiduringpregnancy, patnershivstatus,
         durationofbfmonths_cat, bwt_cat, sexofthebaby, tmembraner_cat) %>%
  tbl_uvregression(
    method       = coxph,
    y            = Surv(surv_time, event),
    exponentiate = TRUE
  ) %>%
  bold_p(t = 0.05)

tbl_cox_merged <- tbl_merge(
  tbls        = list(tbl_cox_mv_sub, tbl_cox_final),
  tab_spanner = c("**Univariable HR**", "**Multivariable HR**")
) %>%
  modify_caption("**Table 11: Cox PH Regression – Combined Table**")

print(tbl_cox_merged)

#spatial analysis


################################################################################
#  SPATIAL MAPPING – HIV-EXPOSED INFANTS, HOMA BAY COUNTY
#  Produces:
#    1. Choropleth map of HIV positivity rate per sub-county
#    2. Dot-density map of individual patient GPS points
#    3. Kernel density (heatmap) of HIV-positive cases
#    4. Faceted maps by adherence, CD4, breastfeeding duration
#    5. Interactive leaflet map (HTML output)

################################################################################

# ── 2. SUB-COUNTY CENTROIDS (approximate, from data medians) ─────────────────
subcounty_coords <- tribble(
  ~subcounty,        ~lon,     ~lat,
  "Homa Bay",       34.704,  -0.495,
  "Kabondo Kasipul",34.475,  -0.527,
  "Karachauony",    34.596,  -0.427,
  "Kasipul",        34.481,  -0.628,
  "Mbita",          34.210,  -0.461,
  "Ndhiwa",         34.508,  -0.538,
  "Ragwe",          34.535,  -0.506,
  "Suba",           34.207,  -0.539
)

# Sub-county summary statistics
sc_summary <- df %>%
  group_by(subcounty) %>%
  summarise(
    n          = n(),
    hiv_pos    = sum(hiv_positive),
    hiv_rate   = mean(hiv_positive) * 100,
    median_age = median(age, na.rm = TRUE),
    good_adherence_pct = mean(adherence == "good", na.rm = TRUE) * 100,
    .groups    = "drop"
  ) %>%
  left_join(subcounty_coords, by = "subcounty")

cat("Sub-county HIV rates:\n")
print(sc_summary %>% select(subcounty, n, hiv_pos, hiv_rate) %>%
        arrange(desc(hiv_rate)))

# ── 3. CONVERT POINTS TO SF OBJECTS ──────────────────────────────────────────
pts_sf <- st_as_sf(df, coords = c("var52", "var51"), crs = 4326)
sc_sf  <- st_as_sf(sc_summary, coords = c("lon", "lat"), crs = 4326)


kenya_counties <- st_read("D:/Desktop/Kenya_Wards/kenya_wards.shp")

plot(kenya_counties)

#filtering homa bay

homa_bay <- kenya_counties %>%
  filter(county == "Homa Bay")


ggplot(homa_bay) +
  geom_sf(fill = "lightblue", color = "black") +
  labs(title = "Homa Bay County") +
  theme_minimal()


library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggspatial)

# ── Ensure CRS matches (VERY IMPORTANT) ──────────────────────
homa_bay <- st_transform(homa_bay, crs = 4326)

# ── 5A. BUBBLE MAP ───────────────────────────────────────────
p_bubble <- ggplot() +
  
  # Background: Homa Bay County
  geom_sf(data = homa_bay, fill = "#f0f0f0",
          color = "black", linewidth = 0.5) +
  
  # All patient points (light grey)
  geom_point(data = df,
             aes(x = var52, y = var51),
             color = "#b0b8c1", size = 0.8, alpha = 0.4) +
  
  # HIV-positive patients (red)
  geom_point(data = df %>% filter(hiv_positive == 1),
             aes(x = var52, y = var51),
             color = "#E63946", size = 2.5, alpha = 0.85) +
  
  # Bubble per sub-county
  geom_point(data = sc_summary,
             aes(x = lon, y = lat, size = n, fill = hiv_rate),
             shape = 21, color = "white", stroke = 1.2, alpha = 0.85) +
  
  # Labels
  geom_text_repel(data = sc_summary,
                  aes(x = lon, y = lat,
                      label = paste0(subcounty, "\n",
                                     round(hiv_rate, 1), "%")),
                  size = 3, fontface = "bold",
                  color = "#1d3557") +
  
  # Scales
  scale_fill_viridis_c(name = "HIV Positivity (%)",
                       option = "plasma", direction = -1) +
  
  scale_size_continuous(name = "Patients (n)",
                        range = c(4, 15)) +
  
  # Map extent (optional but OK)
  coord_sf(xlim = c(33.9, 35.0),
           ylim = c(-1.1, 0.0)) +
  
  # Map elements
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tr",
                         which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  
  # Labels
  labs(
    title = "HIV Positivity Among Exposed Infants — Homa Bay County",
    subtitle = "Bubble size = sample size; colour = HIV positivity rate (%)",
    caption = "Data: Homa Bay PMTCT Programme",
    x = "Longitude", y = "Latitude"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ── Plot ─────────────────────────────────────────────────────
print(p_bubble)



library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(KernSmooth)
library(sf)

# ── Ensure CRS is correct ─────────────────────────────
homa_bay <- st_transform(homa_bay, crs = 4326)

# ── Filter HIV-positive cases ─────────────────────────
hiv_pts <- df %>%
  filter(hiv_positive == 1,
         !is.na(var52),
         !is.na(var51))

# ── Kernel density estimation ─────────────────────────
bw <- c(bandwidth.nrd(hiv_pts$var52),
        bandwidth.nrd(hiv_pts$var51))

dens <- bkde2D(
  cbind(hiv_pts$longitude, hiv_pts$latitude),
  bandwidth = bw * 1.5,
  gridsize = c(150, 150)
)

# Convert to dataframe
dens_df <- expand.grid(lon = dens$x1, lat = dens$x2) %>%
  mutate(density = as.vector(dens$fhat))

# ── Plot ──────────────────────────────────────────────
p_heat <- ggplot() +
  
  # Background: Homa Bay
  geom_sf(data = homa_bay, fill = "#1a1a2e",
          color = "#444", linewidth = 0.4) +
  
  # Density layer
  geom_tile(data = dens_df,
            aes(x = lon, y = lat, fill = density),
            alpha = 0.85) +
  
  # Color scale
  scale_fill_gradientn(
    colors = c("transparent", "#03045e", "#0077b6", "#00b4d8",
               "#90e0ef", "#ffd166", "#ef233c"),
    name = "Case Density"
  ) +
  
  # HIV-positive points
  geom_point(data = hiv_pts,
             aes(x = longitude, y = latitude),
             color = "white", size = 1.5, alpha = 0.9, shape = 3) +
  
  # Subcounty labels
  geom_text_repel(data = sc_summary,
                  aes(x = lon, y = lat, label = subcounty),
                  color = "white", size = 3, fontface = "bold") +
  
  # Map extent
  coord_sf(xlim = c(33.9, 35.0),
           ylim = c(-1.1, 0.0)) +
  
  # Map elements
  annotation_scale(location = "bl", width_hint = 0.25,
                   text_col = "white", line_col = "white") +
  
  annotation_north_arrow(location = "tr",
                         which_north = "true",
                         style = north_arrow_minimal(
                           line_col = "white",
                           text_col = "white"
                         )) +
  
  # Labels
  labs(
    title = "Spatial Density of HIV-Positive Infant Cases",
    subtitle = "Kernel density of confirmed HIV-positive cases",
    caption = "Data: Homa Bay PMTCT Programme",
    x = "Longitude", y = "Latitude"
  ) +
  
  # Theme
  theme_dark(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", color = "white"),
    plot.background = element_rect(fill = "#1a1a2e"),
    panel.background = element_rect(fill = "#1a1a2e")
  )

# ── Plot ──────────────────────────────────────────────
print(p_heat)
