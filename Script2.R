## =========================================================================
## HIV-FREE SURVIVAL AMONG HIV-EXPOSED INFANTS — HOMA BAY COUNTY
## Addresses objectives 1.4.1-1.4.2 and null hypotheses 1.5(1)-1.5(5)
## =========================================================================

## ---- 0. Packages -------------------------------------------------------
pkgs <- c("readxl", "dplyr", "tidyr", "survival", "survminer",
          "broom", "gtsummary", "spdep", "sf", "car")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

## ---- 1. Import & clean --------------------------------------------------
raw <- readxl::read_excel("~/Research_Proposals_Work/Folder/Augustine Gatimu/HOMABAY DATASET.csv.xlsx", sheet = "HOMABAY DATASET")


df <- raw %>%
  mutate(across(where(is.character), ~ trimws(tolower(.x)))) %>%
  mutate(
    ## survival object already Stata-stset: _t = time, _d = event (NA = censored)
    time  = as.numeric(`_t`),
    event = ifelse(is.na(`_d`), 0, `_d`),
    event = as.numeric(event),
    
    ## factors, reference levels chosen a priori
    maritalstatus       = factor(maritalstatus, levels = c("married", "single", "widowed")),
    residence            = factor(residence, levels = c("rural", "urban")),
    employmentstatus     = factor(employmentstatus),
    ancattendance        = factor(ancattendance, levels = c("no", "yes")),
    facilitylevelanc     = factor(facilitylevelanc),
    distancetofacility   = factor(distancetofacility, levels = c("below 1 km", "1-5 km", "above 5 km")),
    hivstatusbeforepregnancy = factor(hivstatusbeforepregnancy),
    syphillis            = factor(syphillis, levels = c("negative", "positive")),
    uti                  = factor(uti, levels = c("no", "yes")),
    miscarriage          = factor(miscarriage, levels = c("no", "yes")),
    abortion             = factor(abortion, levels = c("no", "yes")),
    stillbirth           = factor(stillbirth, levels = c("no", "yes")),
    premature            = factor(premature, levels = c("no", "yes")),
    haartduringpregnancy = factor(haartduringpregnancy, levels = c("no", "yes")),
    cd4cellcountcellsmm3 = factor(cd4cellcountcellsmm3,
                                  levels = c("below 200", "200-500", ">300", "above 500")),
    historyofstiduringpregnancy = factor(historyofstiduringpregnancy, levels = c("no", "yes")),
    genitalwart          = factor(genitalwart, levels = c("no", "yes")),
    vaginaldischargesyndrome = factor(vaginaldischargesyndrome, levels = c("no", "yes")),
    genitalherpes        = factor(genitalherpes, levels = c("no", "yes")),
    treatedduringpregnancy = factor(treatedduringpregnancy, levels = c("no", "yes")),
    malaria              = factor(malaria, levels = c("no", "yes")),
    anaemia              = factor(anaemia, levels = c("no", "yes")),
    hypertention         = factor(hypertention, levels = c("no", "yes")),
    patnershivstatus     = factor(patnershivstatus),
    whohivdiseasestage   = factor(whohivdiseasestage, levels = c("stage i", "stage ii", "stage iii", "stage iv")),
    tmembraner           = factor(tmembraner),
    educationlevel        = factor(educationlevel,
                                   levels = c("never attended school", "primary", "secondary", "college")),
    adherence             = factor(adherence, levels = c("poor", "fair", "good")),
    sexofthebaby          = factor(sexofthebaby, levels = c("female", "male")),
    subcounty             = factor(subcounty),
    
    ## numeric covariates
    age = as.numeric(age), weightkg = as.numeric(weightkg),
    bwtkgs = as.numeric(bwtkgs), babywt6wks = as.numeric(babywt6wks),
    gestationatfirstancwks = as.numeric(gestationatfirstancwks),
    numberofchildren = as.numeric(numberofchildren),
    numberofancvisitsmade = as.numeric(numberofancvisitsmade),
    hemoglobinlevelduringpregnancy = as.numeric(hemoglobinlevelduringpregnancy),
    durationofbfmonths = as.numeric(durationofbfmonths),
    
    ## geo-coordinates for spatial analysis
    lon = as.numeric(var51), lat = as.numeric(var52)
  ) %>%
  filter(!is.na(time), time > 0)

## drop implausible coordinates outside Kenya's bounding box before spatial work
df_geo <- df %>% filter(between(lon, 33.5, 42), between(lat, -5, 5.5))

## ---- 2. Variable groups --------------------------------------------------
maternal_vars <- c("age","weightkg","maritalstatus","residence","employmentstatus",
                   "ancattendance","facilitylevelanc","distancetofacility",
                   "gestationatfirstancwks","hivstatusbeforepregnancy","syphillis","uti",
                   "miscarriage","abortion","stillbirth","premature","haartduringpregnancy",
                   "cd4cellcountcellsmm3","historyofstiduringpregnancy","genitalwart",
                   "vaginaldischargesyndrome","genitalherpes","treatedduringpregnancy",
                   "malaria","anaemia","hypertention","numberofchildren","patnershivstatus",
                   "whohivdiseasestage","numberofancvisitsmade","hemoglobinlevelduringpregnancy",
                   "tmembraner","educationlevel","adherence")

infant_vars <- c("bwtkgs","babywt6wks","durationofbfmonths","sexofthebaby")

## =========================================================================
## OBJECTIVE 1 — HIV-free survival rate of exposed infants
## =========================================================================
surv_obj <- Surv(time = df$time, event = df$event)

km_overall <- survfit(surv_obj ~ 1, data = df)
summary(km_overall)$table                       # median survival, N events
survival_at_t <- summary(km_overall, times = c(6, 12, 24, 52, 78))
print(survival_at_t)                             # HIV-free survival probability at key ages (weeks)

ggsurvplot(km_overall, data = df, risk.table = TRUE, conf.int = TRUE,
           xlab = "Time (weeks)", ylab = "HIV-free survival probability",
           title = "Kaplan-Meier HIV-free survival: HIV-exposed infants, Homa Bay")

## =========================================================================
## OBJECTIVE 2 — Maternal factors associated with HIV antigen positivity
## H0(1): maternal & infant factors not associated with vertical transmission
## H0(2): HR = 1 for covariate x1 adjusted for others (Wald test per covariate)
## =========================================================================
univ_cox <- function(var, data) {
  f <- as.formula(paste("Surv(time, event) ~", var))
  m <- coxph(f, data = data)
  broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>% mutate(variable = var)
}

maternal_univ <- bind_rows(lapply(maternal_vars, univ_cox, data = df))
print(maternal_univ %>% select(variable, term, estimate, conf.low, conf.high, p.value),
      n = Inf)

## carry forward variables significant at p < 0.20 into the multivariable model
sig_maternal <- maternal_univ %>% filter(p.value < 0.20) %>% pull(variable) %>% unique()
f_maternal <- as.formula(paste("Surv(time, event) ~", paste(sig_maternal, collapse = " + ")))
cox_maternal <- coxph(f_maternal, data = df)
summary(cox_maternal)                            # adjusted HRs -> tests H0(1) & H0(2) for maternal set

## =========================================================================
## OBJECTIVE 3 — Infant factors associated with HIV antigen positivity
## =========================================================================
infant_univ <- bind_rows(lapply(infant_vars, univ_cox, data = df))
print(infant_univ %>% select(variable, term, estimate, conf.low, conf.high, p.value), n = Inf)

sig_infant <- infant_univ %>% filter(p.value < 0.20) %>% pull(variable) %>% unique()
f_infant <- as.formula(paste("Surv(time, event) ~", paste(sig_infant, collapse = " + ")))
cox_infant <- coxph(f_infant, data = df)
summary(cox_infant)

## =========================================================================
## OBJECTIVE 4 — Final model of HIV-free survival (maternal + infant, adjusted)
## H0(3): no time-varying (non-proportional) effect of a binary covariate, theta(t)=0
## =========================================================================
final_vars   <- unique(c(sig_maternal, sig_infant))
f_final      <- as.formula(paste("Surv(time, event) ~", paste(final_vars, collapse = " + ")))
cox_full     <- coxph(f_final, data = df)

## backward elimination on AIC for the parsimonious final model
cox_final <- step(cox_full, direction = "backward")
summary(cox_final)

## proportional-hazards / time-varying effect diagnostic -> tests H0(3)
ph_test <- cox.zph(cox_final)
print(ph_test)
plot(ph_test)

## multicollinearity check
car::vif(cox_final)

## forest plot of final adjusted model
ggforest(cox_final, data = df)

## =========================================================================
## SPATIAL ANALYSIS — H0(4): incidences spatially independent
##                     H0(5): no spatial autocorrelation (Moran's I)
## =========================================================================
geo_pts <- df_geo %>% filter(!is.na(lon), !is.na(lat)) %>%
  select(lon, lat, event)

coords <- as.matrix(geo_pts[, c("lon", "lat")])
knn_nb <- spdep::knn2nb(spdep::knearneigh(coords, k = 6))
lw     <- spdep::nb2listw(knn_nb, style = "W")

moran_test <- spdep::moran.test(geo_pts$event, lw)
print(moran_test)                                # p < 0.05 rejects H0(4)/H0(5)

moran_mc <- spdep::moran.mc(geo_pts$event, lw, nsim = 999)
print(moran_mc)

local_moran <- spdep::localmoran(geo_pts$event, lw)
geo_pts$local_I     <- local_moran[, 1]
geo_pts$local_I_pval <- local_moran[, 5]
head(geo_pts[order(geo_pts$local_I_pval), ])      # potential hotspots/clusters

## =========================================================================
## Summary table for write-up (gtsummary)
## =========================================================================
gtsummary::tbl_regression(cox_final, exponentiate = TRUE)