
library(tidyverse)
library(ISLR)
library(rstatix)


get_summary_stats(data = Wage, age, wage, type = "common")


get_summary_stats(data = Wage, age, wage,
                  type = "mean_sd")
get_summary_stats(data = Wage, age,wage,
                  type = "five_number")



#produce a tidy frequency table

bla <- table(Wage$education, Wage$jobclass)


bla |> prop.table(margin =1)


Wage |> freq_table(education, jobclass)


Wage |> freq_table(jobclass, education)


table(Wage$education, Wage$jobclass) |>
  prop.table(margin = 2)
