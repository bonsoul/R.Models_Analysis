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


get_summary_stats(data = Wage, age, wage,
                  type = 'five_number')


get_summary_stats(data = Wage, age, wage,
                  type = 'median')


Wage |> 
  group_by(health, jobclass) |>
  get_summary_stats(
    wage,age, show = c("mean","sd","median","iqr")
  )


##tidy frequency table

bla <- table(Wage$education,Wage$jobclass)
bla


bla |>prop.table(margin = 1)


Wage |> freq_table(education, jobclass)


table(Wage$education, Wage$jobclass) |>
  prop.table(margin =2)



#Assupmtions


library(ggstatsplot)
