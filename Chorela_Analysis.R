# ==============================================================================
# CHOLERA LINE-LIST ANALYSIS — Kenya, 2008–2023
# Descriptive epidemiology + modelling (mortality risk & incidence forecasting)
# ==============================================================================

library(tidyverse)
library(lubridate)
library(lme4)
library(broom)
library(broom.mixed)
library(forecast)
library(scales)

set.seed(123)

# ------------------------------------------------------------------------------
# 1. LOAD & CLEAN
# ------------------------------------------------------------------------------

raw <- read_csv("cleaned_cholera_data.csv", show_col_types = FALSE)

df <- raw %>%
  rename(
    county       = County,
    sub_county   = `Sub County`,
    ward         = Ward,
    village_town = `Village/town`,
    age          = Age,
    sex          = Sex,
    year_raw     = Year,
    date_onset   = `Dates of Onset`,
    date_seen    = `Date seen at Health Facility`,
    outcome      = `Outcome (A/D)`
  ) %>%
  mutate(
    # standardise text fields: trim, collapse whitespace, title case counties
    county     = str_squish(county) %>% str_to_title(),
    county     = case_when(
      county == "Elgeyo-marakwet" ~ "Elgeyo Marakwet",
      county == "Homabay"         ~ "Homa Bay",
      county == "Kajiando"        ~ "Kajiado",
      county == "Tana River"      ~ "Tana River",
      county == "Trans Nzoia"     ~ "Trans Nzoia",
      TRUE ~ county
    ),
    sex        = str_squish(sex) %>% str_to_upper(),
    sex        = na_if(sex, "UNKNOWN"),
    outcome    = str_squish(outcome) %>% str_to_upper(),
    outcome    = case_when(
      outcome %in% c("A", "ALIVE") ~ "Alive",
      outcome %in% c("D", "DEAD")  ~ "Died",
      TRUE ~ NA_character_          # "Unknown" and blanks -> NA
    ),
    date_onset = suppressWarnings(ymd(date_onset)),
    date_seen  = suppressWarnings(ymd(date_seen)),
    # prefer the onset date; fall back to facility date if onset is missing
    event_date = coalesce(date_onset, date_seen),
    year       = year(event_date),
    month      = floor_date(event_date, "month"),
    epiweek    = epiweek(event_date),
    age        = as.numeric(age),
    age_group  = cut(
      age,
      breaks = c(-Inf, 4, 14, 24, 44, 64, Inf),
      labels = c("0-4", "5-14", "15-24", "25-44", "45-64", "65+")
    ),
    delay_days = as.numeric(date_seen - date_onset)
  ) %>%
  filter(!is.na(event_date), year >= 2008, year <= 2023)

cat("Rows after cleaning:", nrow(df), " | Rows dropped:", nrow(raw) - nrow(df), "\n")

# Sense-check the care-seeking delay (onset -> facility) — flag negative/implausible values
cat("\nCare-seeking delay (days), onset to facility:\n")
print(summary(df$delay_days))
df <- df %>% mutate(delay_days = if_else(delay_days < 0 | delay_days > 60, NA_real_, delay_days))

# ------------------------------------------------------------------------------
# 2. DESCRIPTIVE EPIDEMIOLOGY
# ------------------------------------------------------------------------------

## 2a. National epi curve (monthly)
epi_curve <- df %>%
  count(month, name = "cases")

p_epicurve <- ggplot(epi_curve, aes(month, cases)) +
  geom_col(fill = "#c0392b", width = 25) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "Cholera cases by month, Kenya 2008–2023",
       x = NULL, y = "Reported cases") +
  theme_minimal(base_size = 12)

ggsave("epi_curve.png", p_epicurve, width = 10, height = 5, dpi = 150)

## 2b. Case burden by county (top 15)
county_burden <- df %>%
  count(county, sort = TRUE) %>%
  slice_max(n, n = 15)

p_county <- ggplot(county_burden, aes(reorder(county, n), n)) +
  geom_col(fill = "#2980b9") +
  coord_flip() +
  labs(title = "Top 15 counties by reported cholera cases (2008–2023)",
       x = NULL, y = "Cases") +
  theme_minimal(base_size = 12)

ggsave("county_burden.png", p_county, width = 8, height = 6, dpi = 150)

## 2c. Age-sex distribution
age_sex <- df %>%
  filter(!is.na(age_group), !is.na(sex)) %>%
  count(age_group, sex) %>%
  mutate(n = if_else(sex == "M", -n, n))

p_pyramid <- ggplot(age_sex, aes(age_group, n, fill = sex)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = function(x) comma(abs(x))) +
  labs(title = "Age–sex distribution of cholera cases", x = NULL, y = "Cases",
       fill = "Sex") +
  theme_minimal(base_size = 12)

ggsave("age_sex_pyramid.png", p_pyramid, width = 7, height = 5, dpi = 150)

## 2d. Case fatality ratio (CFR): overall and by year
cfr_overall <- df %>%
  filter(!is.na(outcome)) %>%
  summarise(cases = n(), deaths = sum(outcome == "Died"),
            cfr_pct = round(100 * deaths / cases, 2))
cat("\nOverall CFR (excludes 'Unknown' outcomes):\n"); print(cfr_overall)

cfr_by_year <- df %>%
  filter(!is.na(outcome)) %>%
  group_by(year) %>%
  summarise(cases = n(), deaths = sum(outcome == "Died"),
            cfr_pct = round(100 * deaths / cases, 2))

p_cfr_trend <- ggplot(cfr_by_year, aes(year, cfr_pct)) +
  geom_line(color = "#8e44ad", linewidth = 1) +
  geom_point(color = "#8e44ad") +
  labs(title = "Case fatality ratio by year", x = NULL, y = "CFR (%)") +
  theme_minimal(base_size = 12)

ggsave("cfr_trend.png", p_cfr_trend, width = 8, height = 5, dpi = 150)

## 2e. CFR by county (counties with at least 50 classified outcomes)
cfr_by_county <- df %>%
  filter(!is.na(outcome)) %>%
  group_by(county) %>%
  summarise(cases = n(), deaths = sum(outcome == "Died"),
            cfr_pct = round(100 * deaths / cases, 2)) %>%
  filter(cases >= 50) %>%
  arrange(desc(cfr_pct))

cat("\nCFR by county (n>=50 classified outcomes), highest first:\n")
print(cfr_by_county, n = 20)

# ------------------------------------------------------------------------------
# 3. MODEL 1 — Risk factors for death (mixed-effects logistic regression)
#    Outcome: died vs. survived | Fixed effects: age, sex, year
#    Random effect: county (accounts for county-level clustering/reporting
#    differences without burning a degree of freedom per county)
# ------------------------------------------------------------------------------

model_df <- df %>%
  filter(outcome %in% c("Alive", "Died"), !is.na(age), !is.na(sex), age >= 0, age <= 100) %>%
  mutate(
    died      = if_else(outcome == "Died", 1L, 0L),
    age_c     = as.numeric(scale(age, scale = FALSE)),   # centered, keeps coefficients interpretable per year of age
    year_c    = year - median(year),
    sex       = factor(sex, levels = c("F", "M")),
    county    = factor(county)
  )

cat("\nModelling sample size:", nrow(model_df), "\n")

m_mortality <- glmer(
  died ~ age_c + sex + year_c + (1 | county),
  data    = model_df,
  family  = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

cat("\n--- Mixed-effects logistic regression: mortality ---\n")
print(summary(m_mortality))

# Odds ratios with 95% CI — the interpretable version of the coefficients above
or_table <- tidy(m_mortality, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
  select(term, OR = estimate, conf.low, conf.high, p.value)
cat("\nOdds ratios (exponentiated fixed effects):\n")
print(or_table)

# County-level random effects: which counties have unusually high/low mortality
# risk after accounting for age, sex, and year (i.e. beyond case-mix)
county_re <- ranef(m_mortality)$county %>%
  rownames_to_column("county") %>%
  rename(re_logodds = `(Intercept)`) %>%
  arrange(desc(re_logodds))
cat("\nCounties with highest excess mortality risk (random-effect intercepts):\n")
print(head(county_re, 10))
cat("\nCounties with lowest excess mortality risk:\n")
print(tail(county_re, 10))

# ------------------------------------------------------------------------------
# 4. MODEL 2 — Forecasting monthly case counts (time series)
#    ARIMA on the national monthly incidence series, with seasonal terms
# ------------------------------------------------------------------------------

monthly_ts <- epi_curve %>%
  complete(month = seq(min(month), max(month), by = "month"), fill = list(cases = 0)) %>%
  arrange(month)

ts_cases <- ts(monthly_ts$cases,
               start = c(year(min(monthly_ts$month)), month(min(monthly_ts$month))),
               frequency = 12)

fit_arima <- auto.arima(ts_cases, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
cat("\n--- ARIMA model for monthly cholera cases ---\n")
print(summary(fit_arima))

fc <- forecast(fit_arima, h = 12)

p_forecast <- autoplot(fc) +
  labs(title = "Cholera case forecast, next 12 months",
       x = NULL, y = "Reported cases") +
  theme_minimal(base_size = 12)

ggsave("case_forecast.png", p_forecast, width = 10, height = 5, dpi = 150)

cat("\n12-month forecast (point estimate, 80% and 95% intervals):\n")
print(fc)

# Residual diagnostics — check the model isn't leaving structure on the table
cat("\nLjung-Box test on residuals (p > 0.05 supports adequate fit):\n")
print(checkresiduals(fit_arima, plot = FALSE))

# ------------------------------------------------------------------------------
# 5. SAVE KEY TABLES FOR THE WRITE-UP
# ------------------------------------------------------------------------------

write_csv(cfr_by_year,   "cfr_by_year.csv")
write_csv(cfr_by_county, "cfr_by_county.csv")
write_csv(or_table,      "mortality_odds_ratios.csv")
write_csv(as_tibble(fc), "case_forecast_12mo.csv")

cat("\nDone. Outputs written to working directory:\n",
    " - epi_curve.png, county_burden.png, age_sex_pyramid.png, cfr_trend.png, case_forecast.png\n",
    " - cfr_by_year.csv, cfr_by_county.csv, mortality_odds_ratios.csv, case_forecast_12mo.csv\n")

# ------------------------------------------------------------------------------
# 6. SPATIAL HOTSPOT MAP — county-level case burden and CFR
#
#    Boundary source: mikelmaron/kenya-election-data (public domain, IEBC-derived
#    constituency file, dissolved to 48 county-level features). Swap in KNBS/HDX
#    boundaries if you need an authoritative source for publication.
# ------------------------------------------------------------------------------

library(sf)

# Download once and cache locally (re-run is instant after the first time)
geo_path <- "kenya_counties.geojson"
if (!file.exists(geo_path)) {
  download.file(
    "https://raw.githubusercontent.com/mikelmaron/kenya-election-data/master/data/counties.geojson",
    geo_path, quiet = TRUE
  )
}
kenya_sf <- st_read(geo_path, quiet = TRUE)

# Build a normalised join key (strip spaces/hyphens/apostrophes, uppercase) since
# the boundary file and the line list spell a few counties differently
norm_name <- function(x) x %>% str_to_upper() %>% str_remove_all("[^A-Z]")

kenya_sf <- kenya_sf %>%
  mutate(
    county_clean = case_when(
      COUNTY_NAM == "ELEGEYO-MARAKWET" ~ "Elgeyo Marakwet",   # source has a typo
      COUNTY_NAM == "THARAKA - NITHI"  ~ "Tharaka Nithi",
      TRUE ~ str_to_title(COUNTY_NAM)
    ),
    join_key = norm_name(county_clean)
  )

county_stats <- df %>%
  group_by(county) %>%
  summarise(
    cases  = n(),
    deaths = sum(outcome == "Died", na.rm = TRUE),
    classified = sum(!is.na(outcome)),
    cfr_pct = round(100 * deaths / classified, 1)
  ) %>%
  mutate(join_key = norm_name(county))

map_data <- kenya_sf %>%
  left_join(county_stats %>% st_drop_geometry(), by = "join_key")

# Flag any county in the line list that failed to match a boundary (should be
# empty once name variants above are handled)
unmatched <- county_stats %>% anti_join(st_drop_geometry(kenya_sf), by = "join_key")
if (nrow(unmatched) > 0) {
  cat("\nWARNING — counties in the case data with no boundary match:\n")
  print(unmatched$county)
}

## 6a. Case burden choropleth
p_map_cases <- ggplot(map_data) +
  geom_sf(aes(fill = cases), color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "inferno", direction = -1, na.value = "grey90",
                       labels = comma, name = "Cases") +
  labs(title = "Cholera case burden by county, Kenya (2008\u20132023)",
       subtitle = "Grey = no reported cases in this line list") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("hotspot_map_cases.png", p_map_cases, width = 8, height = 8, dpi = 150)

## 6b. CFR choropleth (restricted to counties with a reasonable denominator,
##     so a county with 2 cases and 1 death doesn't paint as an extreme hotspot)
map_data_cfr <- map_data %>%
  mutate(cfr_display = if_else(classified >= 30, cfr_pct, NA_real_))

p_map_cfr <- ggplot(map_data_cfr) +
  geom_sf(aes(fill = cfr_display), color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey90",
                       name = "CFR (%)") +
  labs(title = "Cholera case fatality ratio by county, Kenya (2008\u20132023)",
       subtitle = "Grey = fewer than 30 classified outcomes (CFR unreliable)") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("hotspot_map_cfr.png", p_map_cfr, width = 8, height = 8, dpi = 150)

cat("\nSpatial maps written: hotspot_map_cases.png, hotspot_map_cfr.png\n")