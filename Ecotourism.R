library(dplyr)
library(lubridate)
library(tidytuesdayR)
library(tidyverse)


#loading the datasets

tuesdata <- tidytuesdayR::tt_load(2026, week = 30)

occurrences <- tuesdata$occurrences
tourism <- tuesdata$tourism
weather <- tuesdata$weather

