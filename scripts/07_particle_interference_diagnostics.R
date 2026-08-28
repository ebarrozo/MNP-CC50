#!/usr/bin/env Rscript

# ==============================================================================
# 07_particle_interference_diagnostics.R
#
# PURPOSE
#   Characterize particle-associated fluorescence in the no-cell MNP-only
#   control experiment and test whether it systematically tracks the ORIGINAL
#   EthD-1 apparent cell-death signal.
#
# IMPORTANT
#   The no-cell controls were acquired without the wash steps used before
#   EthD-1 acquisition in cell-containing wells. They therefore characterize
#   potential optical interference but are not a direct estimate of residual
#   post-wash particle signal. No particle-only signal is subtracted from the
#   primary cell-death measurements.
#
# INPUTS
#   data/no_cell/08.26.26-EthD1_mnponly.xlsx
#   results/v5_2/MNP_v5_2_NormalizedExperimentalWells_ALL_QC.csv
#
# OUTPUTS
#   results/interference/v6_5/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# Run from repository root.
project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(project_dir, "scripts")) ||
    !dir.exists(file.path(project_dir, "data"))) {
  stop("Run this script from the repository root.")
}

v5_dir <- file.path(project_dir, "results", "v5_2")
out_dir <- file.path(project_dir, "results", "interference", "v6_5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

well_file <- file.path(v5_dir, "MNP_v5_2_NormalizedExperimentalWells_ALL_QC.csv")
no_cell_file <- file.path(project_dir, "data", "no_cell", "08.26.26-EthD1_mnponly.xlsx")

for (f in c(well_file, no_cell_file)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

write_out <- function(x, filename) {
  readr::write_csv(x, file.path(out_dir, filename))
}

# ------------------------------------------------------------------------------
# 1. Parse the no-cell plate
# ------------------------------------------------------------------------------

nc_raw <- readxl::read_excel(
  no_cell_file,
  sheet = "End point",
  col_names = FALSE,
  .name_repair = "minimal"
)

# readxl can trim leading blank rows. Locate the plate header dynamically.
header_candidates <- which(vapply(seq_len(nrow(nc_raw)), function(i) {
  vals <- suppressWarnings(as.numeric(unlist(nc_raw[i, 2:13], use.names = FALSE)))
  length(vals) == 12 &&
    all(is.finite(vals)) &&
    identical(as.integer(vals), 1:12)
}, logical(1)))

if (length(header_candidates) != 1) {
  stop(
    "Could not uniquely locate the no-cell plate header containing columns 1-12. ",
    "Candidates: ", paste(header_candidates, collapse = ", ")
  )
}

plate_header_row <- header_candidates[1]
plate_data_rows <- (plate_header_row + 1):(plate_header_row + 8)
plate_rows <- toupper(trimws(
  as.character(unlist(nc_raw[plate_data_rows, 1], use.names = FALSE))
))

if (!identical(plate_rows, LETTERS[1:8])) {
  stop("No-cell plate rows beneath the detected header are not A-H.")
}

row_key <- tribble(
  ~plate_row, ~microplastic, ~size,
  "A", "PS",  "Small",
  "B", "PE",  "Small",
  "C", "PVC", "Small",
  "D", "PS",  "Large",
  "E", "PE",  "Large",
  "F", "PVC", "Large",
  "G", "Mix", "Small",
  "H", "Mix", "Large"
)

dose_key <- tibble(
  plate_col = 1:9,
  dose = rep(c(0.2, 2, 20), each = 3),
  no_cell_technical_replicate = rep(1:3, 3)
)

nc_samples <- map_dfr(seq_along(plate_rows), function(i) {
  tibble(
    plate_row = plate_rows[i],
    plate_col = 1:9,
    signal_no_cell = as.numeric(
      unlist(nc_raw[plate_data_rows[i], 2:10], use.names = FALSE)
    )
  )
}) %>%
  left_join(row_key, by = "plate_row") %>%
  left_join(dose_key, by = "plate_col")

if (nrow(nc_samples) != 72 || anyNA(nc_samples$signal_no_cell)) {
  stop("Expected 72 finite no-cell MNP wells (8 particle rows x 9 dose wells).")
}

get_control_triplicate <- function(data_row_index) {
  as.numeric(unlist(nc_raw[data_row_index, 11:13], use.names = FALSE))
}

nc_pbs      <- get_control_triplicate(plate_data_rows[1])
nc_methanol <- get_control_triplicate(plate_data_rows[2])
nc_media    <- get_control_triplicate(plate_data_rows[3])

nc_control_summary <- tibble(
  control = c("PBS", "Methanol", "Media"),
  mean_signal = c(mean(nc_pbs), mean(nc_methanol), mean(nc_media)),
  sd_signal = c(sd(nc_pbs), sd(nc_methanol), sd(nc_media)),
  n = 3L
)

no_cell_media_mean <- mean(nc_media)

nc_condition <- nc_samples %>%
  group_by(microplastic, size, dose) %>%
  summarise(
    n_no_cell = n(),
    no_cell_mean = mean(signal_no_cell),
    no_cell_sd = sd(signal_no_cell),
    no_cell_min = min(signal_no_cell),
    no_cell_max = max(signal_no_cell),
    .groups = "drop"
  ) %>%
  mutate(
    no_cell_media_mean = no_cell_media_mean,
    particle_excess_signed = no_cell_mean - no_cell_media_mean,
    fold_vs_no_cell_media = no_cell_mean / no_cell_media_mean
  ) %>%
  arrange(microplastic, size, dose)

stopifnot(nrow(nc_condition) == 24)

write_out(nc_control_summary, "V6_5_NoCell_ControlSummary.csv")
write_out(nc_samples, "V6_5_NoCell_RawMNPWells.csv")
write_out(nc_condition, "V6_5_NoCell_ParticleFluorescence_ByCondition.csv")

# ------------------------------------------------------------------------------
# 2. Load the original QC-valid cell-containing assay data
# ------------------------------------------------------------------------------

wells <- readr::read_csv(well_file, show_col_types = FALSE)

required <- c(
  "cell_line", "replicate", "microplastic", "size", "dose",
  "timepoint_hr", "Primary_Media_Valid", "Percent_Dead_Media"
)
missing <- setdiff(required, names(wells))
if (length(missing) > 0) {
  stop("Missing required v5.2 columns: ", paste(missing, collapse = ", "))
}

wells_primary <- wells %>%
  mutate(
    cell_line = recode(
      as.character(cell_line),
      "THP1" = "THP-1",
      .default = as.character(cell_line)
    ),
    microplastic = as.character(microplastic),
    size = str_to_title(as.character(size)),
    dose = as.numeric(dose),
    timepoint_hr = as.character(timepoint_hr),
    replicate = as.character(replicate),
    apparent_death = as.numeric(Percent_Dead_Media)
  ) %>%
  filter(Primary_Media_Valid, is.finite(apparent_death))

dose_mean <- wells_primary %>%
  group_by(cell_line, replicate, microplastic, size, dose, timepoint_hr) %>%
  summarise(
    apparent_death_mean = mean(apparent_death),
    apparent_death_sd = if_else(n() >= 2, sd(apparent_death), NA_real_),
    n_technical = n(),
    .groups = "drop"
  )

joined <- dose_mean %>%
  left_join(
    nc_condition %>%
      select(
        microplastic, size, dose, no_cell_mean, no_cell_sd,
        no_cell_media_mean, particle_excess_signed, fold_vs_no_cell_media
      ),
    by = c("microplastic", "size", "dose")
  )

if (anyNA(joined$particle_excess_signed)) {
  stop("At least one cell-containing exposure condition did not match the no-cell plate.")
}

write_out(joined, "V6_5_DoseMean_Joined_NoCellSignal.csv")

# ------------------------------------------------------------------------------
# 3. Condition-level summaries and Spearman diagnostics
# ------------------------------------------------------------------------------

condition_overall <- joined %>%
  group_by(
    microplastic, size, dose,
    particle_excess_signed, no_cell_mean, fold_vs_no_cell_media
  ) %>%
  summarise(
    mean_apparent_death = mean(apparent_death_mean, na.rm = TRUE),
    median_apparent_death = median(apparent_death_mean, na.rm = TRUE),
    n_bio_dose_means = n(),
    .groups = "drop"
  )

write_out(condition_overall, "V6_5_ExposureCondition_OverallSummary.csv")

spearman_test <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

spearman_overall <- bind_rows(
  spearman_test(
    condition_overall$particle_excess_signed,
    condition_overall$mean_apparent_death
  ) %>% mutate(outcome = "Mean apparent percent death", .before = 1),
  spearman_test(
    condition_overall$particle_excess_signed,
    condition_overall$median_apparent_death
  ) %>% mutate(outcome = "Median apparent percent death", .before = 1)
) %>%
  mutate(
    n_conditions = nrow(condition_overall),
    q_BH = p.adjust(p_value, "BH")
  ) %>%
  select(outcome, n_conditions, spearman_rho, p_value, q_BH)

write_out(spearman_overall, "V6_5_Spearman_Overall.csv")

cell_condition <- joined %>%
  group_by(cell_line, microplastic, size, dose, particle_excess_signed) %>%
  summarise(mean_apparent_death = mean(apparent_death_mean), .groups = "drop")

spearman_cell <- cell_condition %>%
  group_by(cell_line) %>%
  group_modify(~ spearman_test(.x$particle_excess_signed, .x$mean_apparent_death) %>%
                 mutate(n_conditions = nrow(.x), .before = 1)) %>%
  ungroup() %>%
  mutate(q_BH = p.adjust(p_value, "BH"))

write_out(spearman_cell, "V6_5_Spearman_ByCellLine.csv")

time_condition <- joined %>%
  group_by(timepoint_hr, microplastic, size, dose, particle_excess_signed) %>%
  summarise(mean_apparent_death = mean(apparent_death_mean), .groups = "drop")

spearman_time <- time_condition %>%
  group_by(timepoint_hr) %>%
  group_modify(~ spearman_test(.x$particle_excess_signed, .x$mean_apparent_death) %>%
                 mutate(n_conditions = nrow(.x), .before = 1)) %>%
  ungroup() %>%
  mutate(q_BH = p.adjust(p_value, "BH"))

write_out(spearman_time, "V6_5_Spearman_ByTimepoint.csv")

dose_condition <- condition_overall %>%
  group_by(dose) %>%
  group_modify(~ spearman_test(.x$particle_excess_signed, .x$mean_apparent_death) %>%
                 mutate(n_conditions = nrow(.x), .before = 1)) %>%
  ungroup() %>%
  mutate(q_BH = p.adjust(p_value, "BH"))

write_out(dose_condition, "V6_5_Spearman_WithinDose.csv")

cell_time_condition <- joined %>%
  group_by(
    cell_line, timepoint_hr, microplastic, size, dose, particle_excess_signed
  ) %>%
  summarise(mean_apparent_death = mean(apparent_death_mean), .groups = "drop")

within_cell_time <- cell_time_condition %>%
  group_by(cell_line, timepoint_hr) %>%
  group_modify(~ spearman_test(.x$particle_excess_signed, .x$mean_apparent_death) %>%
                 mutate(n_conditions = nrow(.x), .before = 1)) %>%
  ungroup() %>%
  mutate(q_BH = p.adjust(p_value, "BH"))

write_out(within_cell_time, "V6_5_WithinCellTime_Correlations.csv")

cat("\n================ PARTICLE-INTERFERENCE DIAGNOSTIC ================\n")
cat("Unique polymer x size x dose conditions: 24\n")
cat("No-cell controls were measured without the cell-assay wash steps.\n\n")
cat("Overall condition-level Spearman correlations:\n")
print(spearman_overall)
cat("\nBy cell line:\n")
print(spearman_cell)
cat("\nBy exposure duration:\n")
print(spearman_time)
cat("\nWithin dose:\n")
print(dose_condition)
cat("==================================================================\n")
