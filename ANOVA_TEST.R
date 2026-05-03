
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


#anova model
model <- aov(recovery ~ drug * dosage, data = dat)
anova_tbl <- tidy(model)



# Display nicely
anova_tbl %>%
  kable(
    caption = "Two-Way ANOVA Results",
    digits = 4,
    col.names = c("Term", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width = FALSE,
    position = "center"
  ) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#2C3E50") %>%
  column_spec(1, bold = TRUE)




library(flextable)

drug_tukey <- TukeyHSD(model, which = "drug")$drug %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Comparison")


flextable(drug_tukey) %>%
  set_caption("Tukey HSD Post Hoc Test: Drug") %>%
  autofit() %>%
  theme_vanilla() %>%
  colformat_double(digits = 4)


# Tukey test for Interaction
interaction_tukey <- TukeyHSD(model, which = "drug:dosage")$`drug:dosage` %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Comparison")

flextable(interaction_tukey) %>%
  set_caption("Tukey HSD Post Hoc Test: Drug × Dosage") %>%
  autofit() %>%
  theme_vanilla() %>%
  colformat_double(digits = 4)

library(effectsize)


# Interaction plot — the crossing/diverging lines show the interaction
interaction.plot(
  x.factor     = dat$dosage,
  trace.factor = dat$drug,
  response     = dat$recovery,
  fun          = mean,
  type         = "b",
  col          = c("steelblue", "tomato", "forestgreen"),
  pch          = c(16, 17, 18),
  lwd          = 2,
  main         = "Interaction Plot: Drug × Dosage",
  xlab         = "Dosage",
  ylab         = "Mean Recovery Score",
  trace.label  = "Drug"
)
