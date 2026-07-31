library(dplyr)
library(lubridate)
library(tidytuesdayR)
library(tidyverse)
library(knitr)
library(kableExtra)


#loading the datasets

tuesdata <- tidytuesdayR::tt_load(2026, week = 30)

occurrences <- tuesdata$occurrences
tourism <- tuesdata$tourism
weather <- tuesdata$weather


occurrences <- occurrences |>
  mutate(date = as.Date(date), quarter = quarter(date))

weather <- weather |>
  mutate(date = as.Date(date), quarter = quarter(date))


##QN 1 : Under which weather conditions aremost likely to observe a Gouldian finch?


finch <- occurrences |>
  filter(organism_name == "Gouldian finch")

finch_weather <- finch |>
  inner_join(weather, by = c("ws_id","date"),suffix = c("_obs","_wx"))


finch_stations <- unique(finch_weather$ws_id)
baseline_weather <- weather |> filter(ws_id %in% finch_stations)




finch_compare <- bind_rows(
  finch_weather |> mutate(group = "Finch sighting day"),
  baseline_weather |> mutate(group = "All days (same stations)")
)

finch_summary <- finch_compare |> 
  group_by(group) |>
  summarise(
    n_days = n(),
    mean_temp = mean(temp, na.rm = TRUE),
    mean_min_temp = mean(min, na.rm = TRUE),
    mean_max_temp = mean(max, na.rm = TRUE),
    mean_rh = mean(rh, na.rm = TRUE),
    mean_prcp = mean(prcp, na.rm = TRUE),
    pct_rainy  =mean(rainy, na.rm = TRUE),
    mean_wind = mean(wind_speed, na.rm = TRUE),
    .groups = "drop"
  )


print(finch_summary)


finch_summary |>
  kable(
    caption = "Summary Statistics for Gouldian Finch Observations",
    digits = 2
  ) |>
  kable_styling(
    bootstrap_options = c("striped", "hover"),
    full_width = FALSE
  ) |>
  row_spec(0, bold = TRUE, color = "white", background = "#2C3E50") |>
  column_spec(1, bold = TRUE)

##visualization
ggplot(finch_compare, aes(x = group, y = temp, fill = group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(
    title = "Temperature: Gouldian finch sighting days vs. baseline",
    subtitle = "Baseline = all days recorded at stations where finches were seen",
    x = NULL, y = "Mean daily temperature (°C)"
  ) +
  theme_minimal()

##visual two

ggplot(finch_compare, aes(x = group, y = rh, fill = group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(
    title = "Relative humidity: Gouldian finch sighting days vs. baseline",
    x = NULL, y = "Relative humidity (%)"
  ) +
  theme_minimal()

##Q2: How does weather affect tourism numbers in each region?

weather_quarterly <- weather |>
  group_by(ws_id,year,quarter) |>
  summarise(
    mean_temp = mean(temp, na.rm = TRUE),
    total_prcp = sum(prcp, na.rm = TRUE),
    pct_rainy_days = mean(rainy, na.rm = TRUE) * 100,
    mean_wind = mean(wind_speed, na.rm=TRUE),
    .groups = "drop"
  )


tourism_region_q <- tourism |>
  filter(!is.na(ws_id)) |>
  group_by(region,ws_id,year,quarter) |>
  summarise(total_trips = sum(trips, na.rm = TRUE), .groups = "drop") |>
  left_join(weather_quarterly, by = c("ws_id", "year","quarter")) |>
  filter(!is.na(mean_temp))


reqion_weather_cor <- tourism_region_q |>
  group_by(region) |>
  filter(n() >= 8) |>
  summarise(
    n_quarters = n(),
    cor_temp = cor(mean_temp, total_trips, use = "complete.obs"),
    cor_prcp = cor(total_prcp, total_trips, use = "complete.obs"),
    cor_rainy_days = cor(pct_rainy_days, total_trips, use = "complete.obs"),
    .groups = "drop"
  ) |>
  arrange(desc(abs(cor_temp)))

print(reqion_weather_cor, n = 25)

##correlation

reqion_weather_cor |>
  slice_head(n = 25) |>
  kable(
    caption = "Weather Correlations by Region (First 25 Rows)",
    digits = 3,
    align = "c"
  ) |>
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "center"
  )

## plot the region with the strongerst temp-trips
top_region <- reqion_weather_cor$region[1]


tourism_region_q |>
  filter(region == top_region) |>
  ggplot(aes(x = mean_temp, y = total_trips)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = paste("Tourism vs Temperature:", top_region),
    x = "Mean quarterly temp (celcius)", y = "Total trips"
  ) +
  theme_minimal()

##regions where rains seems most -ve

reqion_weather_cor |>
  arrange(cor_rainy_days) |>
  slice_head(n = 10) |>
  kable(
    caption = "Top 10 Regions with the Lowest Correlation for Rainy Days",
    digits = 3
  ) |>
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE
  )
