#!/usr/bin/env Rscript

# ==============================================================================
# 08_particle_interference_condition_models.R
#
# PURPOSE
#   Final condition-level diagnostic testing whether particle-only fluorescence
#   is associated with original apparent EthD-1 cell death after adjustment for
#   polymer treatment, nominal particle-size class, and concentration.
#
#   The analysis uses the 24 unique polymer x size x dose conditions. These are
#   unique analytical conditions from the no-cell control experiment and should
#   not be interpreted as 24 independently repeated no-cell experiments.
#
# INPUT
#   results/interference/v6_5/V6_5_ExposureCondition_OverallSummary.csv
#
# OUTPUT
#   results/interference/v6_6/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(project_dir, "scripts")) ||
    !dir.exists(file.path(project_dir, "results"))) {
  stop("Run this script from the repository root.")
}

in_dir <- file.path(project_dir, "results", "interference", "v6_5")
out_dir <- file.path(project_dir, "results", "interference", "v6_6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

infile <- file.path(in_dir, "V6_5_ExposureCondition_OverallSummary.csv")
if (!file.exists(infile)) stop("Missing input: ", infile)

dat <- readr::read_csv(infile, show_col_types = FALSE) %>%
  mutate(
    microplastic = factor(microplastic),
    size = factor(size),
    dose_f = factor(dose, levels = c(0.2, 2, 20)),
    z_particle_excess = as.numeric(scale(particle_excess_signed))
  )

if (nrow(dat) != 24) {
  stop("Expected 24 unique exposure conditions; found ", nrow(dat))
}

fit_mean_base <- lm(
  mean_apparent_death ~ microplastic + size + dose_f,
  data = dat
)
fit_median_base <- lm(
  median_apparent_death ~ microplastic + size + dose_f,
  data = dat
)
fit_mean_particle <- lm(
  mean_apparent_death ~ z_particle_excess + microplastic + size + dose_f,
  data = dat
)
fit_median_particle <- lm(
  median_apparent_death ~ z_particle_excess + microplastic + size + dose_f,
  data = dat
)

mean_coef <- broom::tidy(fit_mean_particle, conf.int = TRUE) %>%
  mutate(outcome = "Mean apparent percent death") %>%
  select(outcome, everything())

median_coef <- broom::tidy(fit_median_particle, conf.int = TRUE) %>%
  mutate(outcome = "Median apparent percent death") %>%
  select(outcome, everything())

readr::write_csv(
  mean_coef,
  file.path(out_dir, "V6_6_ConditionLevel_Model_Mean.csv")
)
readr::write_csv(
  median_coef,
  file.path(out_dir, "V6_6_ConditionLevel_Model_Median.csv")
)

particle_effects <- bind_rows(mean_coef, median_coef) %>%
  filter(term == "z_particle_excess") %>%
  transmute(
    outcome,
    n_conditions = nrow(dat),
    estimate_percent_death_per_1SD_particle_excess = estimate,
    std_error = std.error,
    statistic,
    p_value = p.value,
    ci_lower = conf.low,
    ci_upper = conf.high
  ) %>%
  mutate(q_BH = p.adjust(p_value, "BH"))

readr::write_csv(
  particle_effects,
  file.path(out_dir, "V6_6_ParticleSignal_EffectSummary.csv")
)

cmp_mean <- anova(fit_mean_base, fit_mean_particle) %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  as_tibble() %>%
  mutate(outcome = "Mean apparent percent death")

cmp_median <- anova(fit_median_base, fit_median_particle) %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  as_tibble() %>%
  mutate(outcome = "Median apparent percent death")

readr::write_csv(
  bind_rows(cmp_mean, cmp_median),
  file.path(out_dir, "V6_6_ModelComparison.csv")
)

diag <- tibble(
  outcome = rep(
    c("Mean apparent percent death", "Median apparent percent death"),
    each = 2
  ),
  model = rep(
    c(
      "Polymer + size + dose",
      "Polymer + size + dose + particle-only signal"
    ),
    2
  ),
  n_conditions = 24,
  AIC = c(
    AIC(fit_mean_base), AIC(fit_mean_particle),
    AIC(fit_median_base), AIC(fit_median_particle)
  ),
  BIC = c(
    BIC(fit_mean_base), BIC(fit_mean_particle),
    BIC(fit_median_base), BIC(fit_median_particle)
  ),
  adjusted_R2 = c(
    summary(fit_mean_base)$adj.r.squared,
    summary(fit_mean_particle)$adj.r.squared,
    summary(fit_median_base)$adj.r.squared,
    summary(fit_median_particle)$adj.r.squared
  )
)

readr::write_csv(
  diag,
  file.path(out_dir, "V6_6_ModelDiagnostics.csv")
)

# HC3 robust-SE sensitivity if optional packages are installed.
if (
  requireNamespace("sandwich", quietly = TRUE) &&
  requireNamespace("lmtest", quietly = TRUE)
) {
  robust_one <- function(fit, outcome_label) {
    ct <- lmtest::coeftest(
      fit,
      vcov. = sandwich::vcovHC(fit, type = "HC3")
    )
    tibble(
      outcome = outcome_label,
      term = rownames(ct),
      estimate = ct[, 1],
      robust_SE_HC3 = ct[, 2],
      statistic = ct[, 3],
      p_value = ct[, 4]
    ) %>%
      filter(term == "z_particle_excess")
  }

  robust_effects <- bind_rows(
    robust_one(fit_mean_particle, "Mean apparent percent death"),
    robust_one(fit_median_particle, "Median apparent percent death")
  ) %>%
    mutate(q_BH = p.adjust(p_value, "BH"))

  readr::write_csv(
    robust_effects,
    file.path(out_dir, "V6_6_ParticleSignal_HC3_Sensitivity.csv")
  )
}

cat("\n================ CONDITION-LEVEL INTERFERENCE MODEL ================\n")
cat("Unique exposure conditions: ", nrow(dat), "\n\n", sep = "")
print(particle_effects)
cat("\nModel diagnostics:\n")
print(diag)
cat("====================================================================\n")
