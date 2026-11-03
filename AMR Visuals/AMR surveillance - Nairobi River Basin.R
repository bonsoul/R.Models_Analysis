# ============================================================
# AMR surveillance - Nairobi River Basin (E. coli / Klebsiella)
# Input : amr_merged.csv
# Output: figures/*.png and tables/*.csv
# ============================================================

library(tidyverse)
library(scales)

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# ---------- 0. SETTINGS YOU MAY NEED TO EDIT ----------------
data_file <- "C:/Users/ADMIN/Documents/R.Models_Analysis/amr_outputs/amr_merged.csv"

# The file has no date column, only "round". Map each round to a year here.
# (Both set to 2025 - change if Round 2 was collected in another year.)
round_year <- c("Round 1" = 2025, "Round 2" = 2025)
report_year <- 2025

# Colours
pal_sir <- c(S = "#2E8B57", I = "#F2B134", R = "#C0392B")
pal_org <- c("E. coli"                  = "#1F77B4",
             "Klebsiella pneumoniae/spp." = "#E67E22",
             "K. oxytoca"               = "#7F7F7F")

theme_set(theme_minimal(base_size = 12) +
            theme(plot.title = element_text(face = "bold"),
                  panel.grid.minor = element_blank()))

save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path("figures", paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
}

# ---------- 1. LOAD & CLEAN ---------------------------------
raw <- read_csv(data_file, show_col_types = FALSE)

amr <- raw %>%
  mutate(
    across(where(is.character), str_squish),
    year = unname(round_year[round]),
    # standardise organism names (isolate_raw has typos / mucoid variants)
    organism = case_when(
      str_detect(str_to_lower(pathogen), "no growth")                    ~ "No growth",
      str_detect(str_to_lower(paste(pathogen, isolate_raw)), "oxytoca")  ~ "K. oxytoca",
      str_detect(str_to_lower(paste(pathogen, isolate_raw)), "coli")     ~ "E. coli",
      str_detect(str_to_lower(pathogen), "klebsiella|pneumoniae")        ~ "Klebsiella pneumoniae/spp.",
      TRUE                                                               ~ pathogen
    ),
    specimen = factor(specimen, levels = c("River water", "Influent", "Effluent")),
    # river name extracted from the site description (for the basin plot)
    river = case_when(
      str_detect(site, regex("slaughter|farmers choice|neema|QMP", ignore_case = TRUE)) ~ "Slaughterhouse discharge",
      str_detect(site, regex("WWTP|waste ?water|effluent", ignore_case = TRUE))        ~ "Dandora WWTP",
      str_detect(site, regex("ngong",    ignore_case = TRUE)) ~ "Ngong River",
      str_detect(site, regex("mbagathi", ignore_case = TRUE)) ~ "Mbagathi River",
      str_detect(site, regex("ruirwaka", ignore_case = TRUE)) ~ "Ruirwaka River",
      str_detect(site, regex("ruiru",    ignore_case = TRUE)) ~ "Ruiru River",
      str_detect(site, regex("thiriri",  ignore_case = TRUE)) ~ "Thiririka River",
      str_detect(site, regex("kamiti",   ignore_case = TRUE)) ~ "Kamiti River",
      str_detect(site, regex("muto mawe", ignore_case = TRUE)) ~ "Muto Mawe River",
      str_detect(site, regex("athi",     ignore_case = TRUE)) ~ "Athi River",
      str_detect(site, regex("nairobi",  ignore_case = TRUE)) ~ "Nairobi River",
      TRUE ~ "Other"
    )
  )

abx_cols <- names(amr)[str_starts(names(amr), "int_")]

# Long format: one row per isolate x antibiotic
long <- amr %>%
  filter(organism != "No growth") %>%
  pivot_longer(all_of(abx_cols), names_to = "antibiotic", values_to = "result",
               names_prefix = "int_") %>%
  filter(!is.na(result)) %>%
  mutate(result = factor(result, levels = c("S", "I", "R")))

isolates <- amr %>% filter(organism != "No growth")

# Quick data checks
cat("Rows:", nrow(amr), "| isolates:", nrow(isolates), "\n")
print(count(amr, organism))
if (!any(amr$organism == "K. oxytoca"))
  message("NOTE: no K. oxytoca isolates in this dataset - it will show as zero in plots 9 & 10.")

# ============================================================
# 3. Priority pathogens by specimen type reported in 2025
# ============================================================
t3 <- isolates %>%
  filter(year == report_year) %>%
  count(specimen, organism, name = "isolates") %>%
  group_by(specimen) %>% mutate(pct = isolates / sum(isolates)) %>% ungroup()
write_csv(t3, "tables/03_pathogen_by_specimen_2025.csv")

p3 <- ggplot(t3, aes(specimen, isolates, fill = organism)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.7) +
  geom_text(aes(label = isolates), position = position_dodge(0.7, preserve = "single"),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = pal_org, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = paste("Priority pathogens by specimen type,", report_year),
       x = NULL, y = "Number of isolates")
save_plot(p3, "03_pathogen_by_specimen")

# ============================================================
# 4 & 5. Resistance against antibiotics: E. coli / Klebsiella
# ============================================================
resistance_plot <- function(org, file_id, title) {
  d <- long %>% filter(organism == org)
  
  tab <- d %>%
    count(antibiotic, result, .drop = FALSE) %>%
    group_by(antibiotic) %>%
    mutate(n_tested = sum(n), pct = n / n_tested) %>% ungroup()
  write_csv(tab, sprintf("tables/%s_resistance_%s.csv", file_id, str_replace_all(org, "\\W+", "_")))
  
  ord <- tab %>% filter(result == "R") %>% arrange(pct) %>% pull(antibiotic)
  tab <- tab %>% mutate(antibiotic = factor(antibiotic, levels = ord))
  
  ggplot(tab, aes(pct, antibiotic, fill = result)) +
    geom_col(width = 0.75) +
    geom_text(data = filter(tab, pct >= 0.06),
              aes(label = percent(pct, accuracy = 1)),
              position = position_stack(vjust = 0.5), colour = "white", size = 3.2) +
    scale_fill_manual(values = pal_sir, name = NULL,
                      labels = c(S = "Susceptible", I = "Intermediate", R = "Resistant")) +
    scale_x_continuous(labels = percent, expand = c(0, 0)) +
    labs(title = title,
         subtitle = paste0("n = ", max(tab$n_tested), " isolates"),
         x = "% of isolates", y = NULL)
}

save_plot(resistance_plot("E. coli", "04", "E. coli: antibiotic susceptibility profile"),
          "04_resistance_ecoli")
save_plot(resistance_plot("Klebsiella pneumoniae/spp.", "05",
                          "Klebsiella: antibiotic susceptibility profile"),
          "05_resistance_klebsiella")

# ============================================================
# 6. Culture workloads
# ============================================================
# One row per sample (a sample can yield >1 isolate)
samples <- amr %>%
  group_by(round, sample_id, specimen) %>%
  summarise(growth = any(organism != "No growth"), .groups = "drop") %>%
  mutate(outcome = if_else(growth, "Growth", "No growth"))

t6 <- samples %>%
  count(round, specimen, outcome) %>%
  group_by(round, specimen) %>%
  mutate(total = sum(n)) %>% ungroup()
write_csv(t6, "tables/06_culture_workload.csv")

p6 <- ggplot(t6, aes(specimen, n, fill = outcome)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            colour = "white", size = 3.5) +
  facet_wrap(~round) +
  scale_fill_manual(values = c(Growth = "#1F77B4", `No growth` = "#9AA5B1"), name = NULL) +
  labs(title = "Culture workload by specimen type and round",
       subtitle = paste0("Total samples cultured: ", nrow(samples)),
       x = NULL, y = "Samples cultured")
save_plot(p6, "06_culture_workload")

# ============================================================
# 7. Critical and high priority pathogens
# ============================================================
# ASSUMPTION (edit if your protocol differs), following WHO BPPL logic:
#   CRITICAL = Enterobacterales resistant to 3rd-gen cephalosporin
#              (ceftazidime or cefepime R) OR carbapenem (meropenem R)
#   HIGH     = not critical, but fluoroquinolone-resistant (ciprofloxacin R)
#   Other    = neither
prio <- isolates %>%
  mutate(
    c3g  = int_Ceftazidime == "R" | int_Cefepime == "R",
    cr   = int_Meropenem == "R",
    fqr  = int_Ciprofloxacin == "R",
    priority = case_when(
      coalesce(cr, FALSE) | coalesce(c3g, FALSE) ~ "Critical",
      coalesce(fqr, FALSE)                       ~ "High",
      TRUE                                       ~ "Not priority"
    ),
    priority = factor(priority, levels = c("Critical", "High", "Not priority"))
  )

t7 <- prio %>%
  count(organism, specimen, priority, name = "isolates") %>%
  group_by(organism, specimen) %>% mutate(pct = isolates / sum(isolates)) %>% ungroup()
write_csv(t7, "tables/07_priority_classification.csv")

p7 <- ggplot(t7, aes(specimen, isolates, fill = priority)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = isolates), position = position_stack(vjust = 0.5),
            colour = "white", size = 3.5) +
  facet_wrap(~organism) +
  scale_fill_manual(values = c(Critical = "#C0392B", High = "#E67E22",
                               `Not priority` = "#9AA5B1"), name = NULL) +
  labs(title = "Critical and high priority resistant isolates",
       subtitle = "Critical: 3GC- or carbapenem-resistant | High: fluoroquinolone-resistant",
       x = NULL, y = "Number of isolates")
save_plot(p7, "07_critical_high_priority", w = 10)

# ============================================================
# 8. Trend in E. coli resistance (Round 1 -> Round 2)
# ============================================================
ec <- long %>% filter(organism == "E. coli") %>%
  group_by(round, antibiotic) %>%
  summarise(n = n(), n_R = sum(result == "R"), pct_R = n_R / n, .groups = "drop")

delta <- ec %>%
  select(round, antibiotic, n, n_R, pct_R) %>%
  pivot_wider(names_from = round, values_from = c(n, n_R, pct_R)) %>%
  mutate(change_pp = (`pct_R_Round 2` - `pct_R_Round 1`) * 100,
         trend = case_when(change_pp < 0 ~ "Decrease",
                           change_pp > 0 ~ "Increase", TRUE ~ "No change"),
         p_fisher = pmap_dbl(
           list(`n_R_Round 1`, `n_Round 1`, `n_R_Round 2`, `n_Round 2`),
           \(r1, n1, r2, n2) fisher.test(matrix(c(r1, n1 - r1, r2, n2 - r2), 2))$p.value)) %>%
  arrange(change_pp)
write_csv(delta, "tables/08_ecoli_resistance_trend.csv")
print(delta %>% select(antibiotic, starts_with("pct_R"), change_pp, trend, p_fisher))

decr <- delta %>% filter(trend == "Decrease") %>% pull(antibiotic)

# one label per Round-2 value (antibiotics with identical % share a label)
lab8 <- ec %>% filter(round == "Round 2") %>%
  group_by(round, pct_R) %>%
  summarise(antibiotic = paste(antibiotic, collapse = " / "), .groups = "drop")

p8 <- ec %>%
  mutate(trend = if_else(antibiotic %in% decr, "Decreased", "Increased / unchanged")) %>%
  ggplot(aes(round, pct_R, group = antibiotic, colour = trend)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  geom_text(data = lab8, aes(round, pct_R, label = antibiotic, group = NULL),
            hjust = -0.1, size = 3, show.legend = FALSE, inherit.aes = FALSE) +
  scale_colour_manual(values = c(Decreased = "#2E8B57", `Increased / unchanged` = "#C0392B"),
                      name = NULL) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  scale_x_discrete(expand = expansion(add = c(0.2, 1.8))) +
  labs(title = "E. coli resistance by antibiotic across sampling rounds",
       subtitle = "Only two time points: read as a change, not a trend",
       x = NULL, y = "% resistant")
save_plot(p8, "08_ecoli_resistance_trend", w = 10)

# ============================================================
# 9. Organism distribution along the Nairobi River Basin (river water)
# ============================================================
org3 <- names(pal_org)

t9 <- isolates %>%
  filter(specimen == "River water", organism %in% org3) %>%
  mutate(organism = factor(organism, levels = org3)) %>%
  count(river, organism, .drop = FALSE, name = "isolates") %>%
  group_by(river) %>% filter(sum(isolates) > 0) %>% ungroup()
write_csv(t9, "tables/09_river_basin_distribution.csv")

river_ord <- t9 %>% group_by(river) %>% summarise(t = sum(isolates)) %>%
  arrange(t) %>% pull(river)

p9 <- t9 %>% mutate(river = factor(river, levels = river_ord)) %>%
  ggplot(aes(isolates, river, fill = organism)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.15)) +
  scale_fill_manual(values = pal_org, name = NULL, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, .08))) +
  labs(title = "Distribution of E. coli and Klebsiella along the Nairobi River Basin",
       subtitle = "River water isolates",
       x = "Number of isolates", y = NULL)
save_plot(p9, "09_river_basin_distribution", w = 10)

# Optional: by site and position along the river (upstream / confluence / downstream)
t9b <- isolates %>%
  filter(specimen == "River water", organism %in% org3) %>%
  count(source, organism, name = "isolates")
write_csv(t9b, "tables/09b_river_position_distribution.csv")

# ============================================================
# 10. Organism distribution in effluent samples
# ============================================================
t10 <- isolates %>%
  filter(specimen == "Effluent", organism %in% org3) %>%
  mutate(organism = factor(organism, levels = org3)) %>%
  count(source, organism, .drop = FALSE, name = "isolates")
write_csv(t10, "tables/10_effluent_distribution.csv")

p10 <- ggplot(t10, aes(source, isolates, fill = organism)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.15), width = 0.8) +
  geom_text(aes(label = isolates), position = position_dodge2(0.8, preserve = "single"),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = pal_org, name = NULL, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Distribution of E. coli and Klebsiella in effluent samples",
       x = "Effluent source", y = "Number of isolates")
save_plot(p10, "10_effluent_distribution")

message("Done. See ./figures and ./tables")