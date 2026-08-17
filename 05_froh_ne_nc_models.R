#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  05 - Relationship between inbreeding (FROH) and population size (Nc, Ne)
#
#  Companion code for:
#    Bourbon et al. (2026) Inbreeding and Demographic History of Caribou
#    (Rangifer tarandus) in Western Canada Inferred From Genome-Wide SNP Data.
#    Evolutionary Applications. doi:10.1111/eva.70311
#
#  Pipeline (see Methods 2.6):
#    Model log(FROH) as a function of log(Nc) and log(Ne). Candidate models
#    include a null, linear Nc, linear Ne, additive, quadratic (Nc^2, Ne^2),
#    and an Nc x Ne interaction. Models are ranked by AICc; those within
#    delta AICc < 2 (excluding the null) are conditionally averaged (MuMIn).
#    For the top model(s), the percent change in FROH per 10-fold change in the
#    predictor is reported as (10^beta - 1) * 100%.
#
#    Run at two scales:
#      1. individual subpopulations (n > 5, Nc > 5)
#      2. pooled genetic units + subpopulations not assigned to a unit
#
#  Requires: MuMIn, dplyr, readxl
#  Inputs  : FROH_by_subpopulation.txt              (mean FROH per subpop, script 02)
#            census_size.xlsx                        (Nc per subpopulation)
#            ne_subpopulations.xlsx                  (currentNe2 subpop Ne, script 04)
#            ne_genetic_units.xlsx                   (currentNe2 unit Ne, script 04)
#            sample_sizes.txt                        (genomic n per subpop)
#            genetic_unit_map.txt                    (subpop -> genetic unit)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(MuMIn)
library(dplyr)
library(readxl)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
data_dir <- "data"
out_dir  <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ---------------------------------------------------------------------------
# 1. Load inputs
# ---------------------------------------------------------------------------
# Mean FROH per subpopulation (from script 02)
froh <- read.table(file.path(out_dir, "FROH_by_subpopulation.txt"),
                   header = TRUE, sep = "\t") %>%
  rename(population = group, FROH = mean_FROH) %>%
  select(population, FROH)

# Genomic sample size per subpopulation. Columns: population, n
sample_sizes <- read.table(file.path(data_dir, "sample_sizes.txt"),
                           header = TRUE, sep = "\t")

# Census size per subpopulation. Columns: population, Nc (mature individuals)
census <- read_excel(file.path(data_dir, "census_size.xlsx")) %>%
  select(population, Nc) %>%
  mutate(Nc = as.numeric(Nc))

# Contemporary Ne per subpopulation at r = 1.0 (panmictic; from script 04).
# Columns: population, Ne_point, Ne_lower90, Ne_upper90
ne_subpop <- read_excel(file.path(data_dir, "ne_subpopulations.xlsx")) %>%
  select(population, Ne = Ne_point) %>%
  mutate(Ne = as.numeric(Ne))

# Contemporary Ne per pooled genetic unit at r = 1.0 (from script 04).
ne_units <- read_excel(file.path(data_dir, "ne_genetic_units.xlsx")) %>%
  select(population = genetic_unit, Ne = Ne_point) %>%
  mutate(Ne = as.numeric(Ne))

# Map of subpopulation -> pooled genetic unit. Columns: population, genetic_unit
unit_map_df <- read.table(file.path(data_dir, "genetic_unit_map.txt"),
                          header = TRUE, sep = "\t")
genetic_unit_map <- setNames(unit_map_df$genetic_unit, unit_map_df$population)


# ---------------------------------------------------------------------------
# 2. Model set, selection, averaging, and effect size
# ---------------------------------------------------------------------------
# Fits the candidate models, ranks by AICc, averages models within delta < 2
# (excluding the null), and reports the percent change in FROH per 10-fold
# change in each retained predictor.
run_models <- function(df, label = "") {

  df <- as.data.frame(df)

  models <- list(
    null        = lm(log_FROH ~ 1,                               data = df, na.action = na.fail),
    Nc          = lm(log_FROH ~ log_Nc,                          data = df, na.action = na.fail),
    Ne          = lm(log_FROH ~ log_Ne,                          data = df, na.action = na.fail),
    Nc_Ne       = lm(log_FROH ~ log_Nc + log_Ne,                 data = df, na.action = na.fail),
    Nc_quad     = lm(log_FROH ~ log_Nc + I(log_Nc^2),            data = df, na.action = na.fail),
    Ne_quad     = lm(log_FROH ~ log_Ne + I(log_Ne^2),            data = df, na.action = na.fail),
    Nc_Ne_inter = lm(log_FROH ~ log_Nc * log_Ne,                 data = df, na.action = na.fail)
  )

  aic <- data.frame(
    model  = names(models),
    AICc   = sapply(models, AICc),
    R2     = sapply(models, function(m) summary(m)$r.squared),
    stringsAsFactors = FALSE
  )
  aic$delta  <- aic$AICc - min(aic$AICc)
  aic$weight <- exp(-0.5 * aic$delta) / sum(exp(-0.5 * aic$delta))
  aic <- aic[order(aic$AICc), ]

  cat("\n===== Model selection:", label, "=====\n")
  print(aic, row.names = FALSE)

  cat("\n-- Nc-only model --\n");  print(summary(models$Nc))
  cat("\n-- Ne-only model --\n");  print(summary(models$Ne))

  # Top models within delta < 2, excluding the null
  top <- aic$model[aic$delta < 2]
  top <- top[top != "null"]

  if (length(top) == 1) {
    best <- models[[top]]
    cat("\n===== Single best model:", top, "=====\n")
    print(summary(best))
    coefs <- coef(best)
  } else {
    avg <- model.avg(model.sel(models[top]), revised.var = TRUE, full = FALSE)
    cat("\n===== Model-averaged (conditional):", label, "=====\n")
    print(summary(avg))
    coefs <- coef(avg)
  }

  # Percent change in FROH per 10-fold change in a predictor: (10^beta - 1)*100
  for (term in c("log_Nc", "log_Ne")) {
    if (term %in% names(coefs)) {
      pct <- (10^coefs[[term]] - 1) * 100
      cat(sprintf("Percent change in FROH per 10-fold change in %s: %.1f%%\n",
                  sub("log_", "", term), pct))
    }
  }

  write.table(aic, file.path(out_dir, paste0("FROH_model_selection_",
              gsub("[^A-Za-z0-9]+", "_", label), ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(models = models, aic = aic, top = top))
}


# ---------------------------------------------------------------------------
# 3. Scale 1: individual subpopulations
# ---------------------------------------------------------------------------
# Restrict to subpopulations with genomic n > 5 and a census estimate (Nc > 5).
df_subpop <- froh %>%
  left_join(ne_subpop,    by = "population") %>%
  left_join(census,       by = "population") %>%
  left_join(sample_sizes, by = "population") %>%
  filter(!is.na(Nc), Nc > 5, !is.na(Ne), !is.na(n), n > 5) %>%
  mutate(log_FROH = log(FROH), log_Nc = log(Nc), log_Ne = log(Ne))

cat("Subpopulations in scale-1 analysis:", nrow(df_subpop), "\n")
res_subpop <- run_models(df_subpop, "subpopulations")


# ---------------------------------------------------------------------------
# 4. Scale 2: pooled genetic units + unassigned subpopulations
# ---------------------------------------------------------------------------
# FROH for a unit is the sample-size-weighted mean of its member subpopulations;
# Nc is the sum of member census sizes; Ne comes from the unit-level currentNe2.
assign_unit <- function(pop) {
  u <- genetic_unit_map[pop]
  ifelse(is.na(u), pop, u)
}

froh_units <- froh %>%
  left_join(sample_sizes, by = "population") %>%
  mutate(unit = assign_unit(population)) %>%
  group_by(unit) %>%
  summarise(FROH = weighted.mean(FROH, w = n, na.rm = TRUE),
            n    = sum(n, na.rm = TRUE),
            .groups = "drop") %>%
  rename(population = unit)

nc_units <- census %>%
  mutate(unit = assign_unit(population)) %>%
  group_by(unit) %>%
  summarise(Nc = sum(Nc, na.rm = TRUE), .groups = "drop") %>%
  rename(population = unit)

# Pooled genetic units
df_units <- froh_units %>%
  left_join(ne_units, by = "population") %>%
  left_join(nc_units, by = "population") %>%
  filter(!is.na(Nc), Nc > 5, !is.na(Ne)) %>%
  mutate(log_FROH = log(FROH), log_Nc = log(Nc), log_Ne = log(Ne))

# Subpopulations not assigned to any unit (analysed on their own)
assigned_pops <- names(genetic_unit_map)
df_standalone <- df_subpop %>%
  filter(!population %in% assigned_pops) %>%
  select(population, FROH, n, Nc, Ne, log_FROH, log_Nc, log_Ne)

df_combined <- bind_rows(df_units, df_standalone)

cat("\nScale-2 analysis: ", nrow(df_units), "genetic units +",
    nrow(df_standalone), "standalone subpopulations =",
    nrow(df_combined), "units\n")
res_combined <- run_models(df_combined, "genetic units and standalone")

cat("\nDone. Model-selection tables written to", out_dir, "\n")
