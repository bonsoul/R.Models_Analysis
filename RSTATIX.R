library(tidyverse)
library(ISLR)
library(rstatix)
library(gtsummary)


#descriptive stats

get_summary_stats(data = Wage, age, wage, type = 'common')



tbl_summary(
  Wage,
  include = c(age, wage)
)
