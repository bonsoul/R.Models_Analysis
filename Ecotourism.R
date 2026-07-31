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


ggplot(finch_compare, aes(x = group, y = temp, fill = group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(
    title = "Temperature: Gouldian finch sighting days vs. baseline",
    subtitle = "Baseline = all days recorded at stations where finches were seen",
    x = NULL, y = "Mean daily temperature (°C)"
  ) +
  theme_minimal()


