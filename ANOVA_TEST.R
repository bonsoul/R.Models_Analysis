
library(knitr)
library(kableExtra)
library(broom)
library(dplyr)

set.seed(42)

dat <- data.frame(
  recovery = c(
    rnorm(10, mean = 60, sd = 8),   # Drug A, Low dose
    rnorm(10, mean = 75, sd = 8),   # Drug A, High dose
    rnorm(10, mean = 55, sd = 8),   # Drug B, Low dose
    rnorm(10, mean = 58, sd = 8),   # Drug B, High dose  ← dosage barely helps
    rnorm(10, mean = 65, sd = 8),   # Drug C, Low dose
    rnorm(10, mean = 85, sd = 8)    # Drug C, High dose  ← dosage helps a lot
  ),
  drug    = factor(rep(c("A","B","C"), each = 20)),
  dosage  = factor(rep(rep(c("Low","High"), each = 10), times = 3))
)

head(dat, 12)


aggregate(recovery ~ drug +dosage, data = dat, FUN = mean)


# model

model_check <- aov(recovery ~ drug * dosage, data = dat)



shapiro.test(residuals(model_check))


library(car)
leveneTest(recovery ~ drug * dosage, data = dat)


plot(model_check, which = 1)
plot(model_check, which = 2)



model <- aov(recovery ~ drug * dosage, data = dat)
summary(model)
