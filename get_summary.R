
library(tidyverse)
library(ISLR)
library(rstatix)


get_summary_stats(data = Wage, age, wage, type = "common")


get_summary_stats(data = Wage, age, wage,
                  type = "mean_sd")
get_summary_stats(data = Wage, age,wage,
                  type = "five_number")
