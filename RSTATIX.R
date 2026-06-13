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




ggbetweenstats(
  data = airquality,
  x = Month,
  y = Ozone,
  type = "np",
  pairwise.display = 'none',
  results.subtitle =F
)


airquality |>
  group_by(Month) |>
  identify_outliers(Ozone)


Auto |> shapiro_test(acceleration)


Auto |> 
  group_by(cylinders, origin) |>
  shapiro_test(acceleration,horsepower)


Auto |> 
  group_by(cylinders, origin) |>
  shapiro_test(acceleration,horsepower) |>
  p_round(digits = 2) |>
  p_format(accuracy = 0.01)


## comparing means and ranks

## Mann-Whitney test

data("ToothGrowth")

df <- ToothGrowth |>
  convert_as_factor(dose, supp) |>
  reorder_levels('dose', order = c("2", "1","0.5"))


wilcox_test(
  data = df,
  formula = len ~ supp
)


df |>
  group_by(dose) |>
  wilcox_test(len ~ supp, detailed = TRUE) |>
  adjust_pvalue(method = "bonferroni") |>
  add_significance()
