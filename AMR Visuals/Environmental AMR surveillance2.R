# =============================================================================
# Environmental AMR surveillance, Phase II (2025)
# Sheets used: "ENV DATA - Q2-2025" (Round 1) and "ROUND 2"
# ("Sheet1" is only the antibiotic code/disc-content reference table.)
#
# Outputs (written to ./amr_outputs): one PNG per analysis + CSV tables
# =============================================================================

library(tidyverse)
library(readxl)
library(scales)

path <- "C:/Users/ADMIN/Downloads/ZDU_VEES DOCUMENTS/DATA/ENVIRONMENTAL AMR DATA-phase II -2025.xlsx"
out_dir <- "amr_outputs"
dir.create(out_dir, showWarnings = FALSE)

# ---- 0. Settings ------------------------------------------------------------
abx <- c("Ceftazidime", "Meropenem", "Cefepime", "Cefoxitin",
         "Amoxicillin-clavulanate (AMC)", "Gentamicin", "Ciprofloxacin",
         "Streptomycin", "Cotrimoxazole (SXT)", "Chloramphenicol", "Tetracycline")

# Define colours for SIR categories
pal_sir <- c(
  S = "#2E8B57",  # Susceptible
  I = "#F2B134",  # Intermediate
  R = "#C0392B"   # Resistant
)

# Set the default plotting theme
theme_set(
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
)

# Display the colour palette
pal_sir



save_plot <- function(p, name, w = 9, h = 5.5)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")

# ---- 1. Read the two data sheets --------------------------------------------
# Layout: 11 antibiotics, each with a "Zone" (mm) column followed by an
# interpretation column (R/S/I). Column positions differ slightly by sheet.
read_round <- function(sheet, skip, round_label, cols, use_type) {
  raw <- suppressMessages(read_excel(path, sheet = sheet, skip = skip,
                                     col_names = FALSE, col_types = "text"))
  
  out <- tibble(
    round       = round_label,
    sample_id   = raw[[cols$id]],
    sample_type = if (use_type) raw[[cols$type]] else NA_character_,
    source      = if (!is.na(cols$source)) raw[[cols$source]] else NA_character_,
    site        = raw[[cols$site]],
    isolate_raw = raw[[cols$isolate]],
    synergy     = raw[[cols$synergy]]
  )
  for (k in seq_along(abx)) {
    zc <- cols$first_zone + 2 * (k - 1)
    out[[paste0("zone_", abx[k])]] <- suppressWarnings(as.numeric(raw[[zc]]))
    out[[paste0("int_",  abx[k])]] <- str_to_upper(str_trim(raw[[zc + 1]]))
  }
  out
}

r1 <- read_round("ENV DATA - Q2-2025", skip = 10, round_label = "Round 1", use_type = TRUE,
                 cols = list(id = 1, type = 8, source = 9, site = 10,
                             isolate = 16, first_zone = 17, synergy = 39))
r2 <- read_round("ROUND 2", skip = 9, round_label = "Round 2", use_type = FALSE,
                 cols = list(id = 1, type = NA, source = 8, site = 9,
                             isolate = 14, first_zone = 15, synergy = 37))
# (Round 2 has no usable Sample Type column, so it is derived from the site name below.)

# ---- 2. Clean ---------------------------------------------------------------
amr <- bind_rows(r1, r2) %>%
  filter(!(is.na(sample_id) & is.na(site) & is.na(isolate_raw))) %>%   # blank spacer rows
  group_by(round) %>%
  # a 2nd isolate from the same sample sits on its own row with blank sample fields
  fill(sample_id, sample_type, source, site, .direction = "down") %>%
  ungroup() %>%
  mutate(
    iso_l = str_to_lower(isolate_raw),
    pathogen = case_when(
      str_detect(iso_l, "no growth")  ~ "No growth",
      str_detect(iso_l, "coli")       ~ "E. coli",
      str_detect(iso_l, "oxytoca")    ~ "K. oxytoca",
      str_detect(iso_l, "pneumoniae") ~ "K. pneumoniae",
      str_detect(iso_l, "klebsiella") ~ "Klebsiella spp.",   # not speciated on the form
      TRUE                            ~ NA_character_),
    mucoid = str_detect(iso_l, "mucoid"),
    organism = case_when(pathogen == "E. coli" ~ "E. coli",
                         str_detect(pathogen, "Klebsiella|K\\.") ~ "Klebsiella",
                         TRUE ~ NA_character_),
    site_l = str_to_lower(site),
    # Specimen type: given in Round 1; derived from the site name in Round 2
    specimen = case_when(
      !is.na(sample_type)                      ~ str_to_sentence(str_trim(sample_type)),
      str_detect(site_l, "river|swamp")        ~ "River water",
      str_detect(site_l, "influent|entry")     ~ "Influent",
      TRUE                                     ~ "Effluent"),
    facility = case_when(
      specimen == "River water"                              ~ NA_character_,
      str_detect(site_l, "dandora|wwtp|effluent discharge")  ~ "Dandora WWTP",
      str_detect(site_l, "farmers")                          ~ "Farmers Choice",
      str_detect(site_l, "neema")                            ~ "Neema Livestock",
      str_detect(site_l, "qmp")                              ~ "QMP slaughterhouse",
      str_detect(site_l, "ndumb")                            ~ "Ndumbuini slaughterhouse",
      str_detect(site_l, "kiamaiko")                         ~ "Kiamaiko slaughterhouse",
      TRUE                                                   ~ "Other"),
    river = case_when(
      specimen != "River water"                  ~ NA_character_,
      str_detect(site_l, "ngong")                ~ "Ngong River",
      str_detect(site_l, "mbagathi")             ~ "Mbagathi River",
      str_detect(site_l, "thiriri")              ~ "Thiririka River",
      str_detect(site_l, "ruirwaka|ruiru")       ~ "Ruiru/Ruirwaka River",
      str_detect(site_l, "kamiti")               ~ "Kamiti River",
      str_detect(site_l, "muto mawe")            ~ "Muto Mawe River",
      str_detect(site_l, "\\bathi river")        ~ "Athi River",
      str_detect(site_l, "nairobi")              ~ "Nairobi River",
      TRUE                                       ~ "Other")
  ) %>%
  filter(!is.na(pathogen))   # drops "header" rows whose isolate was written on the next line

# Isolates that were recorded but have no AST results at all (check against raw data!)
no_ast <- amr %>%
  filter(pathogen != "No growth", if_all(starts_with("int_"), is.na)) %>%
  select(round, sample_id, site, pathogen)
if (nrow(no_ast) > 0) {
  message("NOTE: isolates recorded without any AST results (kept in counts, excluded from % resistant):")
  print(no_ast)
}

isolates <- amr %>% filter(pathogen != "No growth")

# Long format for AST results (one row per isolate x antibiotic)
ast_long <- isolates %>%
  select(round, sample_id, specimen, river, facility, site, pathogen, organism,
         synergy, starts_with("int_")) %>%
  pivot_longer(starts_with("int_"), names_to = "antibiotic", names_prefix = "int_",
               values_to = "result") %>%
  filter(result %in% c("S", "I", "R")) %>%
  mutate(antibiotic = factor(antibiotic, levels = abx))

# =============================================================================
# 3. Priority pathogens by specimen type
# =============================================================================
pathogen_levels <- c("E. coli", "K. pneumoniae", "Klebsiella spp.", "K. oxytoca")
spec_levels     <- c("River water", "Influent", "Effluent")

t3 <- isolates %>%
  mutate(pathogen = factor(pathogen, pathogen_levels),
         specimen = factor(specimen, spec_levels)) %>%
  count(specimen, pathogen, .drop = FALSE) %>%
  group_by(specimen) %>% mutate(pct_of_specimen = round(100 * n / sum(n), 1)) %>% ungroup()
write_csv(t3, file.path(out_dir, "03_pathogen_by_specimen.csv"))

p3 <- ggplot(t3, aes(specimen, n, fill = pathogen)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  geom_text(aes(label = n), position = position_dodge2(width = 0.9, preserve = "single"),
            vjust = -0.4, size = 3.5) +
  scale_fill_brewer(palette = "Set2", drop = FALSE) +
  labs(title = "Priority pathogens by specimen type, 2025",
       x = NULL, y = "Number of isolates", fill = "Pathogen")
save_plot(p3, "03_pathogen_by_specimen")

# =============================================================================
# 4 & 5. Resistance against antibiotics: E. coli (4) and Klebsiella (5)
# =============================================================================
resistance_table <- function(org) {
  ast_long %>% filter(organism == org) %>%
    group_by(antibiotic) %>%
    summarise(n_tested = n(), n_R = sum(result == "R"), n_I = sum(result == "I"),
              n_S = sum(result == "S"), pct_R = round(100 * n_R / n_tested, 1),
              .groups = "drop") %>%
    rowwise() %>%
    mutate(ci_low  = round(100 * binom.test(n_R, n_tested)$conf.int[1], 1),
           ci_high = round(100 * binom.test(n_R, n_tested)$conf.int[2], 1)) %>%
    ungroup() %>% arrange(desc(pct_R))
}

plot_resistance <- function(org, tag) {
  tab <- resistance_table(org)
  write_csv(tab, file.path(out_dir, paste0(tag, "_resistance_", str_replace_all(str_to_lower(org), "[^a-z]", ""), ".csv")))
  n_iso <- max(tab$n_tested)   # isolates with AST results
  
  d <- ast_long %>% filter(organism == org) %>%
    count(antibiotic, result) %>%
    group_by(antibiotic) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
    mutate(antibiotic = factor(antibiotic, levels = rev(tab$antibiotic)),
           result = factor(result, c("S", "I", "R")))
  
  p <- ggplot(d, aes(pct, antibiotic, fill = result)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = ifelse(pct >= 7, paste0(round(pct), "%"), "")),
              position = position_stack(vjust = 0.5), colour = "white", size = 3.3) +
    scale_fill_manual(values = pal_sir, name = "Interpretation",
                      labels = c(S = "Susceptible", I = "Intermediate", R = "Resistant")) +
    scale_x_continuous(labels = label_percent(scale = 1), expand = c(0, 0)) +
    labs(title = paste("Antibiotic resistance:", org),
         subtitle = paste0("Both rounds combined (", n_iso, " isolates)"),
         x = "% of isolates", y = NULL)
  save_plot(p, paste0(tag, "_resistance_", str_replace_all(str_to_lower(org), "[^a-z]", "")))
  invisible(tab)
}
tab4 <- plot_resistance("E. coli",    "04")
tab5 <- plot_resistance("Klebsiella", "05")

# =============================================================================
# 6. Culture workload
# =============================================================================
samples <- amr %>%
  group_by(round, sample_id, specimen) %>%
  summarise(n_isolates = sum(pathogen != "No growth"), .groups = "drop") %>%
  mutate(outcome = if_else(n_isolates > 0, "Growth", "No growth"))

t6 <- samples %>%
  group_by(round, specimen) %>%
  summarise(samples_cultured = n(), with_growth = sum(outcome == "Growth"),
            no_growth = sum(outcome == "No growth"), isolates_recovered = sum(n_isolates),
            pct_growth = round(100 * with_growth / samples_cultured, 1), .groups = "drop") %>%
  bind_rows(
    samples %>% group_by(round) %>%
      summarise(specimen = "All specimens", samples_cultured = n(),
                with_growth = sum(outcome == "Growth"), no_growth = sum(outcome == "No growth"),
                isolates_recovered = sum(n_isolates),
                pct_growth = round(100 * with_growth / samples_cultured, 1), .groups = "drop"))
write_csv(t6, file.path(out_dir, "06_culture_workload.csv"))

p6 <- samples %>%
  mutate(specimen = factor(specimen, spec_levels), outcome = factor(outcome, c("Growth", "No growth"))) %>%
  count(round, specimen, outcome) %>%
  ggplot(aes(specimen, n, fill = outcome)) +
  geom_col() +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), colour = "white") +
  facet_wrap(~round) +
  scale_fill_manual(values = c(Growth = "#2C7FB8", `No growth` = "#9AA5B1")) +
  labs(title = "Culture workload: samples cultured and outcome",
       x = NULL, y = "Number of samples", fill = NULL)
save_plot(p6, "06_culture_workload", w = 9, h = 5)

# =============================================================================
# 7. Critical and high priority pathogens (WHO Bacterial Priority Pathogens List)
# For Enterobacterales (E. coli, Klebsiella) the listed resistance phenotypes are:
#   CRITICAL: carbapenem-resistant                      -> meropenem R
#   CRITICAL: 3rd-generation cephalosporin-resistant    -> ceftazidime R
# Fluoroquinolone-resistant Enterobacterales are not on the WHO 2024 list, so no
# "High" category applies to these two organisms. If your protocol defines one,
# add a line to case_when() below, e.g.
#   int_Ciprofloxacin == "R" ~ "High: fluoroquinolone-resistant",
# (checked after the Critical rules, so an isolate is counted once, in its highest tier)
# =============================================================================
prio <- isolates %>%
  mutate(
    priority = case_when(
      int_Meropenem   == "R" ~ "Critical: carbapenem-resistant",
      int_Ceftazidime == "R" ~ "Critical: 3rd-gen cephalosporin-resistant",
      TRUE                   ~ "Not in a priority category"),
    priority = factor(priority, c("Critical: carbapenem-resistant",
                                  "Critical: 3rd-gen cephalosporin-resistant",
                                  "Not in a priority category")))

t7 <- prio %>% count(specimen, organism, priority, .drop = FALSE) %>% arrange(specimen, organism)
write_csv(t7, file.path(out_dir, "07_priority_pathogens.csv"))

p7 <- prio %>%
  mutate(specimen = factor(specimen, spec_levels)) %>%
  count(specimen, organism, priority) %>%
  ggplot(aes(specimen, n, fill = priority)) +
  geom_col() +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), colour = "white", size = 3.3) +
  facet_wrap(~organism) +
  scale_fill_manual(values = c("#7B1E1E", "#D6604D", "#BFC5CC"), drop = FALSE) +
  labs(title = "Critical priority resistant isolates by specimen (WHO BPPL)",
       x = NULL, y = "Number of isolates", fill = NULL) +
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
save_plot(p7, "07_priority_pathogens", w = 10, h = 6)

# =============================================================================
# 8. E. coli resistance trend, Round 1 -> Round 2
# (Only two sampling rounds exist, so the trend is a 2-point comparison.)
# =============================================================================
trend <- ast_long %>% filter(organism == "E. coli") %>%
  group_by(round, antibiotic) %>%
  summarise(n = n(), pct_R = 100 * mean(result == "R"), .groups = "drop")

trend_wide <- trend %>% select(-n) %>%
  pivot_wider(names_from = round, values_from = pct_R) %>%
  mutate(change_pp = `Round 2` - `Round 1`) %>% arrange(change_pp)
write_csv(trend_wide, file.path(out_dir, "08_ecoli_resistance_trend.csv"))

# Multiple antibiotic resistance (MAR) index per isolate = antibiotics R / antibiotics tested
mar <- ast_long %>% filter(organism == "E. coli") %>%
  group_by(round, sample_id, site, pathogen) %>%
  summarise(mar = mean(result == "R"), .groups = "drop")
write_csv(mar %>% group_by(round) %>%
            summarise(n_isolates = n(), mean_MAR = round(mean(mar), 3), median_MAR = round(median(mar), 3)),
          file.path(out_dir, "08_ecoli_MAR_by_round.csv"))

p8 <- ggplot(trend, aes(round, pct_R, group = antibiotic, colour = antibiotic)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 100)) +
  labs(title = "E. coli resistance level by sampling round",
       x = NULL, y = "% resistant", colour = NULL)
save_plot(p8, "08_ecoli_resistance_trend", w = 9.5, h = 5.5)

# =============================================================================
# 9a. Distribution of E. coli / Klebsiella along the Nairobi river basin
# 9b. Distribution in effluent samples
# =============================================================================
org_cols <- c("E. coli" = "#2C7FB8", "K. pneumoniae" = "#D95F0E",
              "Klebsiella spp." = "#F2B134", "K. oxytoca" = "#7A4EAB")

# 9a: river water, per river
t9a <- isolates %>% filter(specimen == "River water") %>%
  mutate(pathogen = factor(pathogen, pathogen_levels)) %>%
  count(river, pathogen, .drop = FALSE)
write_csv(t9a, file.path(out_dir, "09a_river_basin_distribution.csv"))

p9a <- t9a %>%
  mutate(river = fct_reorder(river, n, .fun = sum)) %>%
  ggplot(aes(n, river, fill = pathogen)) +
  geom_col() +
  scale_fill_manual(values = org_cols, drop = FALSE) +
  labs(title = "Microorganisms along the Nairobi river basin",
       subtitle = "River-water isolates, both rounds combined",
       x = "Number of isolates", y = NULL, fill = NULL)
save_plot(p9a, "09a_river_basin_distribution")

# 9a (detail): per sampling site, split by round
p9a2 <- isolates %>% filter(specimen == "River water") %>%
  mutate(pathogen = factor(pathogen, pathogen_levels)) %>%
  count(round, river, site, pathogen) %>%
  ggplot(aes(n, str_wrap(site, 38), fill = pathogen)) +
  geom_col() +
  facet_wrap(~round, scales = "free_y") +
  scale_fill_manual(values = org_cols, drop = FALSE) +
  labs(title = "River-water isolates by sampling site", x = "Number of isolates", y = NULL, fill = NULL)
save_plot(p9a2, "09a_river_sites_by_round", w = 12, h = 8)

# 9b: effluent samples, per facility
t9b <- isolates %>% filter(specimen == "Effluent") %>%
  mutate(pathogen = factor(pathogen, pathogen_levels)) %>%
  count(facility, pathogen, .drop = FALSE)
write_csv(t9b, file.path(out_dir, "09b_effluent_distribution.csv"))

p9b <- t9b %>%
  mutate(facility = fct_reorder(facility, n, .fun = sum)) %>%
  ggplot(aes(n, facility, fill = pathogen)) +
  geom_col() +
  scale_fill_manual(values = org_cols, drop = FALSE) +
  labs(title = "Microorganisms in effluent samples",
       subtitle = "WWTP and slaughterhouse effluent, both rounds combined",
       x = "Number of isolates", y = NULL, fill = NULL)
save_plot(p9b, "09b_effluent_distribution")

message("Done. Tables and figures are in: ", normalizePath(out_dir))

