#!/usr/bin/env Rscript

# =============================================================================
# 01_preprocess_plate_data.R
# Enrico Barrozo, Ph.D.
# Updated 2026-08-26
#
# PURPOSE
# Rebuild the MNP EthD-1 cytotoxicity dataset from the complete set of raw BMG
# MARS workbooks for all five cell lines, using the corrected exposure doses
# (0.2, 2, 20 ug/mL), then calculate replicate-level CC50 values.
#
# EXPECTED INPUT DESIGN
#   5 cell lines x 3 biological replicates x 4 timepoints = 60 raw .xlsx files
#   4 polymers x 2 size classes x 3 doses x 3 technical wells per plate
#
# RAW PLATE MAP (validated against key_plate_layout.xlsx)
#   Rows A-H:
#     A PS Small      B PE Small      C PVC Small
#     D PS Large      E PE Large      F PVC Large
#     G Mix Small     H Mix Large
#   Columns 1-3:  0.2 ug/mL, technical replicates 1-3
#   Columns 4-6:  2   ug/mL, technical replicates 1-3
#   Columns 7-9: 20   ug/mL, technical replicates 1-3
#   Control triplicates:
#     PBS      = row A, columns 10-12
#     Methanol = row B, columns 10-12
#     Media    = row C, columns 10-12
#
# PRIMARY NORMALIZATION
#   QC-valid plates only: Percent_Dead_Media = 100 * (sample - media) / (methanol - media).
#   Failed media-control plates are retained in audit files but primary normalized values are NA.
#   Same-plate PBS and median-control normalizations are generated as sensitivity analyses.
#   No individual well is automatically removed as an outlier.
#
# CC50
#   Technical wells are first averaged within biological replicate x dose.
#   CC50 is then estimated separately for each biological replicate by linear
#   interpolation across adjacent tested doses bracketing 50% death.
#   Boundary values are assay-range/censored values:
#     <=0.2 ug/mL -> stored numerically as 0.2 with censor='Lower'
#     >20 ug/mL   -> stored numerically as 20  with censor='Upper'
#
# No-cell MNP-only fluorescence controls are intentionally NOT incorporated in
# this version; they should be analyzed later as a dedicated interference/QC
# sensitivity analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(stringr)
})

# =============================================================================
# 0. PATHS + SETTINGS
# =============================================================================

# Run this script from the repository root.
project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(project_dir, "scripts")) ||
    !dir.exists(file.path(project_dir, "metadata"))) {
  stop("Run this script from the repository root (the directory containing scripts/ and metadata/).")
}


raw_dir <- file.path(project_dir, "data", "raw")
out_dir <- file.path(project_dir, "results", "v5_2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DOSE_LEVELS <- c(0.2, 2, 20)
TP_LEVELS   <- c(3, 6, 24, 48)
CELL_LEVELS <- c("BeWo", "HTR8", "JEG3", "THP-1", "Jurkat")
POLY_LEVELS <- c("PS", "PE", "PVC", "Mix")
SIZE_LEVELS <- c("Small", "Large")
LOW_CAP     <- min(DOSE_LEVELS)
HIGH_CAP    <- max(DOSE_LEVELS)

# Set TRUE only if you intentionally want output files in the repository root as well.
# Default is FALSE so the raw-input directory remains clean.
WRITE_COMPATIBILITY_COPIES_TO_DATA_DIR <- FALSE

plate_map <- tribble(
  ~plate_row, ~microplastic, ~size,
  "A",        "PS",          "Small",
  "B",        "PE",          "Small",
  "C",        "PVC",         "Small",
  "D",        "PS",          "Large",
  "E",        "PE",          "Large",
  "F",        "PVC",         "Large",
  "G",        "Mix",         "Small",
  "H",        "Mix",         "Large"
)

control_rows <- c(
  PBS      = "A",
  Methanol = "B",
  Media    = "C"
)

# Plate-layout key supplied with the raw data. This is treated as an audit key:
# the key retains the original nominal labels 1/10/100, while the corrected
# actual concentrations used analytically are 0.2/2/20 ug/mL.
plate_key_file <- file.path(project_dir, "metadata", "key_plate_layout.xlsx")

validate_plate_key <- function(path) {
  if (!file.exists(path)) {
    stop("Plate-layout key not found: ", path)
  }

  key <- readxl::read_excel(path, sheet = 1, col_names = FALSE, .name_repair = "minimal")
  if (nrow(key) < 9 || ncol(key) < 13) {
    stop("Plate-layout key does not contain the expected 9 x 13 layout: ", path)
  }

  key <- as.data.frame(key[1:9, 1:13], stringsAsFactors = FALSE)
  expected_rows <- c("A","B","C","D","E","F","G","H")
  observed_rows <- as.character(key[2:9, 1])
  if (!identical(observed_rows, expected_rows)) {
    stop("Plate-layout key row labels do not match A-H. Observed: ",
         paste(observed_rows, collapse = ", "))
  }

  # Validate treatment/size identity from representative wells in column 1.
  expected_prefix <- c("PS-S", "PE-S", "PVC-S", "PS-L",
                       "PE-L", "PVC-L", "Mix-S", "Mix-L")
  observed_col1 <- as.character(key[2:9, 2])
  observed_prefix <- stringr::str_extract(observed_col1, "^[A-Za-z]+-[SL]")
  if (!identical(observed_prefix, expected_prefix)) {
    stop(
      "Plate-layout key treatment/size map does not match the v5 map.\n",
      "Expected: ", paste(expected_prefix, collapse = ", "), "\n",
      "Observed: ", paste(observed_prefix, collapse = ", ")
    )
  }

  # Validate that columns 1-3, 4-6 and 7-9 are triplicate dose blocks.
  for (r in 2:9) {
    x <- as.character(unlist(key[r, 2:10], use.names = FALSE))
    if (!(length(unique(x[1:3])) == 1 &&
          length(unique(x[4:6])) == 1 &&
          length(unique(x[7:9])) == 1)) {
      stop("Plate-layout key does not contain triplicate dose blocks in row ", key[r,1])
    }
  }

  # Validate control positions exactly.
  ctrl <- list(
    PBS      = as.character(unlist(key[2, 11:13], use.names = FALSE)),
    Methanol = as.character(unlist(key[3, 11:13], use.names = FALSE)),
    Media    = as.character(unlist(key[4, 11:13], use.names = FALSE))
  )
  if (!all(ctrl$PBS == "PBS") ||
      !all(ctrl$Methanol == "Methanol") ||
      !all(ctrl$Media == "Media")) {
    stop("Plate-layout key control positions do not match A10:A12 PBS, B10:B12 Methanol, C10:C12 Media.")
  }

  message("Plate-layout key validated: ", basename(path))
  message("  A PS-S; B PE-S; C PVC-S; D PS-L; E PE-L; F PVC-L; G Mix-S; H Mix-L")
  message("  Controls: A10-12 PBS; B10-12 Methanol; C10-12 Media")
  invisible(TRUE)
}

validate_plate_key(plate_key_file)

# =============================================================================
# 1. HELPERS
# =============================================================================

standardize_cell <- function(x) {
  case_when(
    str_detect(x, regex("^bewo$", ignore_case = TRUE)) ~ "BeWo",
    str_detect(x, regex("^htr[-_]?8$", ignore_case = TRUE)) ~ "HTR8",
    str_detect(x, regex("^jeg[-_]?3$", ignore_case = TRUE)) ~ "JEG3",
    str_detect(x, regex("^thp[-_]?1$", ignore_case = TRUE)) ~ "THP-1",
    str_detect(x, regex("^jurkat$", ignore_case = TRUE)) ~ "Jurkat",
    TRUE ~ x
  )
}

parse_mars_filename <- function(path) {
  b <- basename(path)

  # Accepts both:
  # EthD-1 MARS_6hr_btmoptic_Jurkat_3.xlsx
  # EthD-1 MARS_6hr_THP-1_3.xlsx
  m <- str_match(
    b,
    "^EthD-1\\s+MARS_([0-9]+)hr(?:_btmoptic)?_([^_]+)_([123])\\.xlsx$"
  )

  if (all(is.na(m))) {
    stop("Filename did not match expected MARS format: ", b)
  }

  tibble(
    timepoint_hr = as.integer(m[, 2]),
    cell_line    = standardize_cell(m[, 3]),
    replicate    = as.integer(m[, 4]),
    plate_id     = tools::file_path_sans_ext(b),
    source_file  = b
  )
}

read_mars_plate <- function(path) {
  # Raw BMG MARS exports used in this project contain the plate matrix after the
  # row labelled "Raw Data (485, 620)". Reading A10:M18 captures the column
  # labels plus rows A-H for the known exports.
  x <- suppressMessages(
    read_excel(
      path,
      sheet = "End point",
      range = "A10:M18",
      col_names = FALSE
    )
  )

  if (nrow(x) != 9 || ncol(x) != 13) {
    stop(
      "Unexpected raw plate dimensions in ", basename(path),
      ": expected 9 x 13 from A10:M18, observed ", nrow(x), " x ", ncol(x)
    )
  }

  # First row is plate-column labels; remaining 8 rows are A-H.
  p <- x[-1, , drop = FALSE]
  names(p) <- c("plate_row", as.character(1:12))
  p <- p %>% mutate(plate_row = as.character(plate_row))

  if (!identical(p$plate_row, LETTERS[1:8])) {
    stop(
      "Unexpected row labels in ", basename(path), ": ",
      paste(p$plate_row, collapse = ", ")
    )
  }

  p
}

parse_one_plate <- function(path) {
  meta  <- parse_mars_filename(path)
  plate <- read_mars_plate(path)

  # ---------------------------------------------------------------------------
  # Experimental wells: 8 rows x 9 columns = 72 wells per plate
  # ---------------------------------------------------------------------------
  samples <- plate %>%
    select(plate_row, all_of(as.character(1:9))) %>%
    pivot_longer(
      cols = -plate_row,
      names_to = "plate_col",
      values_to = "signal"
    ) %>%
    mutate(
      plate_col = as.integer(plate_col),
      signal = suppressWarnings(as.numeric(signal)),
      dose = DOSE_LEVELS[((plate_col - 1L) %/% 3L) + 1L],
      technical_replicate = ((plate_col - 1L) %% 3L) + 1L
    ) %>%
    left_join(plate_map, by = "plate_row") %>%
    mutate(
      cell_line    = meta$cell_line,
      timepoint_hr = meta$timepoint_hr,
      replicate    = meta$replicate,
      plate_id     = meta$plate_id,
      source_file  = meta$source_file,
      control      = FALSE
    )

  if (anyNA(samples$microplastic) || anyNA(samples$size)) {
    stop("Plate-map join failed for: ", basename(path))
  }

  # ---------------------------------------------------------------------------
  # Controls: rows A/B/C x columns 10:12
  # ---------------------------------------------------------------------------
  controls <- imap_dfr(control_rows, function(row_letter, control_type) {
    vals <- plate %>%
      filter(plate_row == row_letter) %>%
      select(all_of(as.character(10:12))) %>%
      unlist(use.names = FALSE) %>%
      suppressWarnings(as.numeric())

    if (length(vals) != 3 || any(!is.finite(vals))) {
      stop(
        "Could not parse three finite ", control_type,
        " controls from ", basename(path)
      )
    }

    tibble(
      cell_line = meta$cell_line,
      timepoint_hr = meta$timepoint_hr,
      replicate = meta$replicate,
      plate_id = meta$plate_id,
      source_file = meta$source_file,
      control_type = control_type,
      technical_replicate = 1:3,
      control_signal = vals
    )
  })

  list(samples = samples, controls = controls)
}

estimate_cc50 <- function(dose, response) {
  z <- tibble(
    dose = as.numeric(dose),
    response = as.numeric(response)
  ) %>%
    filter(is.finite(dose), is.finite(response)) %>%
    group_by(dose) %>%
    summarise(response = mean(response), .groups = "drop") %>%
    arrange(dose)

  if (nrow(z) < 2) {
    return(tibble(
      cc50 = NA_real_,
      cc50_censor = "Insufficient",
      cc50_label = NA_character_,
      monotonic_nondecreasing = NA,
      crossing_interval = NA_character_
    ))
  }

  monotonic_flag <- all(diff(z$response) >= 0)

  # If 50% death is already reached at the lowest tested dose, true CC50 is
  # at or below the assay lower boundary, regardless of later non-monotonicity.
  if (z$response[1] >= 50) {
    return(tibble(
      cc50 = LOW_CAP,
      cc50_censor = "Lower",
      cc50_label = paste0("<=", LOW_CAP),
      monotonic_nondecreasing = monotonic_flag,
      crossing_interval = paste0("<=", LOW_CAP)
    ))
  }

  # Find the first ADJACENT upward crossing in increasing-dose order.
  # This follows the manuscript definition of interpolation between adjacent
  # concentrations bracketing 50% and remains well-defined for non-monotonic
  # curves such as 20% -> 70% -> 30%.
  up_cross <- which(
    z$response[-nrow(z)] < 50 & z$response[-1] >= 50
  )

  if (length(up_cross) > 0) {
    i <- up_cross[1]
    d1 <- z$dose[i]
    d2 <- z$dose[i + 1]
    y1 <- z$response[i]
    y2 <- z$response[i + 1]

    cc <- if (y2 == y1) {
      mean(c(d1, d2))
    } else {
      d1 + (d2 - d1) * (50 - y1) / (y2 - y1)
    }

    return(tibble(
      cc50 = cc,
      cc50_censor = "None",
      cc50_label = sprintf("%.4g", cc),
      monotonic_nondecreasing = monotonic_flag,
      crossing_interval = paste0(d1, "-", d2)
    ))
  }

  # If no upward crossing occurred and every observed response stayed below 50,
  # the true CC50 is above the tested range.
  if (all(z$response < 50)) {
    return(tibble(
      cc50 = HIGH_CAP,
      cc50_censor = "Upper",
      cc50_label = paste0(">", HIGH_CAP),
      monotonic_nondecreasing = monotonic_flag,
      crossing_interval = paste0(">", HIGH_CAP)
    ))
  }

  # Other non-monotonic patterns that cannot be assigned by the rules above are
  # retained explicitly rather than silently extrapolated.
  tibble(
    cc50 = NA_real_,
    cc50_censor = "Unresolved_nonmonotonic",
    cc50_label = NA_character_,
    monotonic_nondecreasing = monotonic_flag,
    crossing_interval = NA_character_
  )
}

write_both_if_requested <- function(df, filename) {
  write_csv(df, file.path(out_dir, filename))
  if (WRITE_COMPATIBILITY_COPIES_TO_DATA_DIR) {
    write_csv(df, file.path(project_dir, filename))
  }
}

# =============================================================================
# 2. DISCOVER + AUDIT INPUT FILES
# =============================================================================

mars_files <- list.files(
  raw_dir,
  pattern = "^EthD-1\\s+MARS_[0-9]+hr(?:_btmoptic)?_[^_]+_[123]\\.xlsx$",
  full.names = TRUE
)
mars_files <- mars_files[!grepl("^~\\$", basename(mars_files))]

if (length(mars_files) == 0) {
  stop("No raw MARS .xlsx files found in: ", raw_dir)
}

input_manifest <- map_dfr(mars_files, parse_mars_filename) %>%
  arrange(cell_line, replicate, timepoint_hr)

write_csv(input_manifest, file.path(out_dir, "MNP_v5_InputManifest.csv"))

expected_grid <- expand_grid(
  cell_line = CELL_LEVELS,
  replicate = 1:3,
  timepoint_hr = TP_LEVELS
)

missing_inputs <- anti_join(
  expected_grid,
  input_manifest %>% distinct(cell_line, replicate, timepoint_hr),
  by = c("cell_line", "replicate", "timepoint_hr")
)

unexpected_inputs <- anti_join(
  input_manifest %>% distinct(cell_line, replicate, timepoint_hr),
  expected_grid,
  by = c("cell_line", "replicate", "timepoint_hr")
)

duplicate_inputs <- input_manifest %>%
  count(cell_line, replicate, timepoint_hr, name = "n_files") %>%
  filter(n_files != 1)

if (nrow(missing_inputs) > 0) {
  print(missing_inputs)
  stop("Missing expected raw MARS inputs. See table above.")
}
if (nrow(unexpected_inputs) > 0) {
  print(unexpected_inputs)
  stop("Unexpected cell line / replicate / timepoint inputs detected.")
}
if (nrow(duplicate_inputs) > 0) {
  print(duplicate_inputs)
  stop("Duplicate raw plates detected for one or more expected conditions.")
}
if (length(mars_files) != 60) {
  stop("Expected exactly 60 raw MARS workbooks; found ", length(mars_files), ".")
}

message("Input audit passed: 60/60 expected raw MARS plates present.")

# =============================================================================
# 3. PARSE ALL RAW PLATES
# =============================================================================

parsed <- map(mars_files, parse_one_plate)

samples_raw <- map_dfr(parsed, "samples")
controls_raw <- map_dfr(parsed, "controls")

# Expected counts:
# 60 plates * 72 experimental wells = 4,320 rows
# 60 plates * 9 control wells       =   540 rows
stopifnot(nrow(samples_raw) == 60 * 72)
stopifnot(nrow(controls_raw) == 60 * 9)

write_csv(samples_raw, file.path(out_dir, "MNP_v5_RawExperimentalWells_Long.csv"))
write_csv(controls_raw, file.path(out_dir, "MNP_v5_RawControlWells_Long.csv"))

# =============================================================================

# =============================================================================
# 4. CONTROL-WELL QC, ROBUST SUMMARIES, AND PLATE-LEVEL QC CLASSIFICATION
# =============================================================================

# QC thresholds are deliberately configurable and are used to FLAG plates/wells,
# not to delete individual technical wells automatically.
# Primary media-normalized analysis requires only a positive signed dynamic range:
#   Methanol > Media.
# SSMD/Z-prime are QC annotations, not hard primary exclusions. A stricter
# sensitivity analysis additionally requires SSMD_Media >= STRICT_SSMD_MIN.
# This avoids discarding plates solely because of an exploratory separation cutoff
# estimated from only three control wells.
STRICT_SSMD_MIN  <- 1.0
WARN_SSMD_MIN    <- 2.0
TECH_RANGE_WARN_PCTDEAD <- 50
RAW_SIGNAL_REL_RANGE_WARN <- 0.20

# ---- Per-control-well diagnostics ------------------------------------------------
# With n=3 control wells, there is not enough information for reliable automatic
# outlier deletion. We therefore identify the most deviant well from the triplicate
# median and quantify discordance, but retain all wells.
control_well_qc <- controls_raw %>%
  group_by(cell_line, replicate, timepoint_hr, plate_id, source_file, control_type) %>%
  mutate(
    control_mean = mean(control_signal, na.rm = TRUE),
    control_median = median(control_signal, na.rm = TRUE),
    control_sd = sd(control_signal, na.rm = TRUE),
    control_mad = mad(control_signal, center = control_median, constant = 1, na.rm = TRUE),
    control_min = min(control_signal, na.rm = TRUE),
    control_max = max(control_signal, na.rm = TRUE),
    control_range = control_max - control_min,
    abs_dev_from_median = abs(control_signal - control_median),
    max_abs_dev_from_median = max(abs_dev_from_median, na.rm = TRUE),
    most_deviant_well = abs_dev_from_median == max_abs_dev_from_median,
    rel_abs_dev_from_median = if_else(
      is.finite(control_median) & control_median != 0,
      abs_dev_from_median / abs(control_median),
      NA_real_
    )
  ) %>%
  ungroup()

write_csv(control_well_qc, file.path(out_dir, "MNP_v5_2_ControlWell_QC.csv"))

# ---- Plate-level mean and median control summaries -------------------------------
control_summary_long <- controls_raw %>%
  group_by(cell_line, replicate, timepoint_hr, plate_id, source_file, control_type) %>%
  summarise(
    control_mean = mean(control_signal, na.rm = TRUE),
    control_median = median(control_signal, na.rm = TRUE),
    control_sd = sd(control_signal, na.rm = TRUE),
    control_mad = mad(control_signal, center = control_median, constant = 1, na.rm = TRUE),
    control_min = min(control_signal, na.rm = TRUE),
    control_max = max(control_signal, na.rm = TRUE),
    control_range = control_max - control_min,
    n_control = sum(is.finite(control_signal)),
    .groups = "drop"
  )

control_summary <- control_summary_long %>%
  pivot_wider(
    names_from = control_type,
    values_from = c(control_mean, control_median, control_sd, control_mad,
                    control_min, control_max, control_range, n_control),
    names_sep = "_"
  ) %>%
  rename(
    Fmin_PBS = control_mean_PBS,
    Fmax = control_mean_Methanol,
    Fmin_Media = control_mean_Media,
    Fmin_PBS_median = control_median_PBS,
    Fmax_median = control_median_Methanol,
    Fmin_Media_median = control_median_Media,
    SD_PBS = control_sd_PBS,
    SD_Methanol = control_sd_Methanol,
    SD_Media = control_sd_Media,
    MAD_PBS = control_mad_PBS,
    MAD_Methanol = control_mad_Methanol,
    MAD_Media = control_mad_Media,
    Range_PBS = control_range_PBS,
    Range_Methanol = control_range_Methanol,
    Range_Media = control_range_Media,
    N_PBS = n_control_PBS,
    N_Methanol = n_control_Methanol,
    N_Media = n_control_Media
  ) %>%
  mutate(
    DynamicRange_PBS = Fmax - Fmin_PBS,
    DynamicRange_Media = Fmax - Fmin_Media,
    DynamicRange_PBS_median = Fmax_median - Fmin_PBS_median,
    DynamicRange_Media_median = Fmax_median - Fmin_Media_median,

    SSMD_Media = if_else(
      sqrt(SD_Methanol^2 + SD_Media^2) > 0,
      DynamicRange_Media / sqrt(SD_Methanol^2 + SD_Media^2),
      NA_real_
    ),
    SSMD_PBS = if_else(
      sqrt(SD_Methanol^2 + SD_PBS^2) > 0,
      DynamicRange_PBS / sqrt(SD_Methanol^2 + SD_PBS^2),
      NA_real_
    ),
    Zprime_Media = if_else(
      abs(DynamicRange_Media) > 0,
      1 - 3 * (SD_Methanol + SD_Media) / abs(DynamicRange_Media),
      NA_real_
    ),
    Zprime_PBS = if_else(
      abs(DynamicRange_PBS) > 0,
      1 - 3 * (SD_Methanol + SD_PBS) / abs(DynamicRange_PBS),
      NA_real_
    ),

    Valid_DynamicRange_Media = is.finite(DynamicRange_Media) & DynamicRange_Media > 0,
    Valid_DynamicRange_PBS = is.finite(DynamicRange_PBS) & DynamicRange_PBS > 0,

    Media_QC = case_when(
      !Valid_DynamicRange_Media ~ "FAIL_direction",
      !is.finite(SSMD_Media) ~ "WARN_missing_separation",
      SSMD_Media < STRICT_SSMD_MIN ~ "WARN_poor_separation",
      SSMD_Media < WARN_SSMD_MIN ~ "WARN_separation",
      TRUE ~ "PASS"
    ),
    PBS_QC = case_when(
      !Valid_DynamicRange_PBS ~ "FAIL_direction",
      !is.finite(SSMD_PBS) ~ "WARN_missing_separation",
      SSMD_PBS < STRICT_SSMD_MIN ~ "WARN_poor_separation",
      SSMD_PBS < WARN_SSMD_MIN ~ "WARN_separation",
      TRUE ~ "PASS"
    ),

    # Primary validity is based only on correct control direction.
    Primary_Media_Valid = Valid_DynamicRange_Media,
    PBS_Sensitivity_Valid = Valid_DynamicRange_PBS,
    # Stricter sensitivity sets additionally require SSMD >= 1.
    Strict_Media_Valid = Valid_DynamicRange_Media & is.finite(SSMD_Media) & SSMD_Media >= STRICT_SSMD_MIN,
    Strict_PBS_Valid = Valid_DynamicRange_PBS & is.finite(SSMD_PBS) & SSMD_PBS >= STRICT_SSMD_MIN,

    # Same-plate rescue candidates. These are for sensitivity analyses only.
    Median_Media_Rescue_Candidate = !Primary_Media_Valid &
      is.finite(DynamicRange_Media_median) & DynamicRange_Media_median > 0,
    PBS_Rescue_Candidate = !Primary_Media_Valid & PBS_Sensitivity_Valid,

    Plate_QC_Interpretation = case_when(
      Primary_Media_Valid & Media_QC == "PASS" ~ "Primary media valid",
      Primary_Media_Valid & Media_QC == "WARN_separation" ~ "Primary media valid; weak separation",
      Primary_Media_Valid & Media_QC == "WARN_poor_separation" ~ "Primary media valid; poor separation flag",
      Primary_Media_Valid & Media_QC == "WARN_missing_separation" ~ "Primary media valid; separation metric unavailable",
      !Primary_Media_Valid & PBS_Rescue_Candidate & Median_Media_Rescue_Candidate ~ "Primary media failed; PBS and median-media sensitivity available",
      !Primary_Media_Valid & PBS_Rescue_Candidate ~ "Primary media failed; PBS sensitivity available",
      !Primary_Media_Valid & Median_Media_Rescue_Candidate ~ "Primary media failed; median-media sensitivity available",
      TRUE ~ "Primary media failed; no same-plate rescue"
    )
  )

write_csv(control_summary, file.path(out_dir, "MNP_v5_2_PlateControl_QC.csv"))

write_csv(
  control_summary %>%
    count(cell_line, Media_QC, name = "n_plates") %>%
    arrange(cell_line, Media_QC),
  file.path(out_dir, "MNP_v5_2_PlateControl_QC_byCellLine.csv")
)

# =============================================================================
# 5. NORMALIZATION: PRIMARY MEAN-MEDIA + SAME-PLATE SENSITIVITY ROUTES
# =============================================================================

normalized <- samples_raw %>%
  left_join(
    control_summary %>%
      select(
        cell_line, replicate, timepoint_hr, plate_id,
        Fmax, Fmin_PBS, Fmin_Media,
        Fmax_median, Fmin_PBS_median, Fmin_Media_median,
        DynamicRange_PBS, DynamicRange_Media,
        DynamicRange_PBS_median, DynamicRange_Media_median,
        SSMD_Media, SSMD_PBS, Zprime_Media, Zprime_PBS,
        Media_QC, PBS_QC, Primary_Media_Valid, PBS_Sensitivity_Valid, Strict_Media_Valid, Strict_PBS_Valid,
        Median_Media_Rescue_Candidate, PBS_Rescue_Candidate,
        Plate_QC_Interpretation
      ),
    by = c("cell_line", "replicate", "timepoint_hr", "plate_id")
  ) %>%
  mutate(
    # Raw formula outputs are retained for audit regardless of plate QC.
    Percent_Dead_Media_raw_unfiltered = 100 * (signal - Fmin_Media) / DynamicRange_Media,
    Percent_Dead_PBS_raw_unfiltered = 100 * (signal - Fmin_PBS) / DynamicRange_PBS,
    Percent_Dead_MediaMedian_raw_unfiltered = 100 * (signal - Fmin_Media_median) / DynamicRange_Media_median,

    # Primary media values are NA on failed plates. No cross-plate controls are borrowed.
    Percent_Dead_Media_raw = if_else(
      Primary_Media_Valid,
      Percent_Dead_Media_raw_unfiltered,
      NA_real_
    ),
    Percent_Dead_Media = if_else(
      Primary_Media_Valid,
      pmin(pmax(Percent_Dead_Media_raw_unfiltered, 0), 100),
      NA_real_
    ),

    # PBS sensitivity is only considered valid when same-plate PBS QC passes.
    Percent_Dead_PBS_raw = if_else(
      PBS_Sensitivity_Valid,
      Percent_Dead_PBS_raw_unfiltered,
      NA_real_
    ),
    Percent_Dead_PBS = if_else(
      PBS_Sensitivity_Valid,
      pmin(pmax(Percent_Dead_PBS_raw_unfiltered, 0), 100),
      NA_real_
    ),

    # Median-control sensitivity: same plate, same control type, robust center.
    # This never replaces the primary mean-control analysis automatically.
    Percent_Dead_MediaMedian = if_else(
      is.finite(DynamicRange_Media_median) & DynamicRange_Media_median > 0,
      pmin(pmax(Percent_Dead_MediaMedian_raw_unfiltered, 0), 100),
      NA_real_
    ),

    Methanol_PBS = if_else(PBS_Sensitivity_Valid, 100, NA_real_),
    Media_PBS = if_else(PBS_Sensitivity_Valid,
                        100 * (Fmin_Media - Fmin_PBS) / DynamicRange_PBS,
                        NA_real_),
    PBS_PBS = if_else(PBS_Sensitivity_Valid, 0, NA_real_),
    Methanol_Media = if_else(Primary_Media_Valid, 100, NA_real_),
    Media_Media = if_else(Primary_Media_Valid, 0, NA_real_),
    PBS_Media = if_else(Primary_Media_Valid,
                        100 * (Fmin_PBS - Fmin_Media) / DynamicRange_Media,
                        NA_real_)
  )

if (anyNA(normalized$Fmax) || anyNA(normalized$Fmin_Media)) {
  stop("Control join failed for at least one experimental well.")
}

write_csv(normalized, file.path(out_dir, "MNP_v5_2_NormalizedExperimentalWells_ALL_QC.csv"))

# =============================================================================
# 6. SERVER-COMPATIBLE NORMALIZED CSVs (PRIMARY QC-AWARE MEDIA)
# =============================================================================

normalized_compat <- normalized %>%
  transmute(
    signal = signal,
    microplastic = microplastic,
    size = size,
    dose = dose,
    control = FALSE,
    cell_line = cell_line,
    timepoint = paste0(timepoint_hr, "hr"),
    replicate = replicate,
    plate_id = plate_id,
    Fmax = Fmax,
    Fmin_PBS = Fmin_PBS,
    Fmin_Media = Fmin_Media,
    Percent_Dead_PBS = Percent_Dead_PBS,
    Percent_Dead_Media = Percent_Dead_Media,
    Methanol_PBS = Methanol_PBS,
    Media_PBS = Media_PBS,
    PBS_PBS = PBS_PBS,
    Methanol_Media = Methanol_Media,
    Media_Media = Media_Media,
    PBS_Media = PBS_Media,
    Media_QC = Media_QC,
    PBS_QC = PBS_QC,
    Primary_Media_Valid = Primary_Media_Valid,
    PBS_Sensitivity_Valid = PBS_Sensitivity_Valid,
    SSMD_Media = SSMD_Media,
    SSMD_PBS = SSMD_PBS
  ) %>%
  arrange(cell_line, replicate, timepoint, microplastic, size, dose)

write_csv(
  normalized_compat,
  file.path(out_dir, "Rawdata_Merged_QC_PercentDeadPBSandMedia_NormByPlate_ALL_v5_2.csv")
)

for (cl in CELL_LEVELS) {
  x <- normalized_compat %>% filter(cell_line == cl)
  stopifnot(nrow(x) == 12 * 72)
  fn <- paste0("Rawdata_Merged_QC_PercentDeadPBSandMedia_NormByPlate_", cl, "_v5_2.csv")
  write_csv(x, file.path(out_dir, fn))
}

# =============================================================================
# 7. TECHNICAL-REPLICATE QC (FLAG ONLY; NO AUTOMATIC WELL EXCLUSION)
# =============================================================================

tech_qc <- normalized %>%
  group_by(cell_line, replicate, timepoint_hr, plate_id, microplastic, size, dose) %>%
  summarise(
    n_tech = n_distinct(technical_replicate),
    n_finite_primary_media = sum(is.finite(Percent_Dead_Media)),
    raw_signal_mean = mean(signal, na.rm = TRUE),
    raw_signal_median = median(signal, na.rm = TRUE),
    raw_signal_sd = sd(signal, na.rm = TRUE),
    raw_signal_min = min(signal, na.rm = TRUE),
    raw_signal_max = max(signal, na.rm = TRUE),
    raw_signal_range = raw_signal_max - raw_signal_min,
    raw_signal_rel_range = if_else(
      is.finite(raw_signal_median) & raw_signal_median != 0,
      raw_signal_range / abs(raw_signal_median),
      NA_real_
    ),
    mean_dead_media = if_else(
      any(is.finite(Percent_Dead_Media)),
      mean(Percent_Dead_Media, na.rm = TRUE),
      NA_real_
    ),
    median_dead_media = if_else(
      any(is.finite(Percent_Dead_Media)),
      median(Percent_Dead_Media, na.rm = TRUE),
      NA_real_
    ),
    sd_dead_media = if_else(
      sum(is.finite(Percent_Dead_Media)) >= 2,
      sd(Percent_Dead_Media, na.rm = TRUE),
      NA_real_
    ),
    min_dead_media = { x <- Percent_Dead_Media[is.finite(Percent_Dead_Media)]; if (length(x)) min(x) else NA_real_ },
    max_dead_media = { x <- Percent_Dead_Media[is.finite(Percent_Dead_Media)]; if (length(x)) max(x) else NA_real_ },
    range_dead_media = if_else(is.finite(min_dead_media) & is.finite(max_dead_media), max_dead_media - min_dead_media, NA_real_),
    Media_QC = first(Media_QC),
    Primary_Media_Valid = first(Primary_Media_Valid),
    .groups = "drop"
  ) %>%
  mutate(
    Flag_raw_signal_discordance = is.finite(raw_signal_rel_range) & raw_signal_rel_range >= RAW_SIGNAL_REL_RANGE_WARN,
    Flag_pctdead_discordance = is.finite(range_dead_media) & range_dead_media >= TECH_RANGE_WARN_PCTDEAD,
    Technical_QC_Flag = case_when(
      !Primary_Media_Valid ~ "PLATE_QC_FAILED",
      Flag_raw_signal_discordance & Flag_pctdead_discordance ~ "WARN_raw_and_normalized_discordance",
      Flag_raw_signal_discordance ~ "WARN_raw_signal_discordance",
      Flag_pctdead_discordance ~ "WARN_normalized_discordance",
      TRUE ~ "PASS"
    )
  )

write_csv(tech_qc, file.path(out_dir, "MNP_v5_2_TechnicalReplicate_QC.csv"))

# Per-well experimental deviation from triplicate median for manual inspection.
experimental_well_qc <- normalized %>%
  group_by(cell_line, replicate, timepoint_hr, plate_id, microplastic, size, dose) %>%
  mutate(
    triplicate_signal_median = median(signal, na.rm = TRUE),
    abs_dev_signal_from_median = abs(signal - triplicate_signal_median),
    max_abs_dev_signal_from_median = max(abs_dev_signal_from_median, na.rm = TRUE),
    most_deviant_signal_well = abs_dev_signal_from_median == max_abs_dev_signal_from_median
  ) %>%
  ungroup() %>%
  select(
    cell_line, replicate, timepoint_hr, plate_id, microplastic, size, dose,
    technical_replicate, signal, triplicate_signal_median,
    abs_dev_signal_from_median, most_deviant_signal_well,
    Percent_Dead_Media, Media_QC, Primary_Media_Valid
  )

write_csv(experimental_well_qc, file.path(out_dir, "MNP_v5_2_ExperimentalWell_QC.csv"))

# =============================================================================
# 8. BIOLOGICAL-REPLICATE DOSE MEANS: PRIMARY + SENSITIVITY NORMALIZATIONS
# =============================================================================

bio_dose_summary <- normalized %>%
  group_by(cell_line, replicate, microplastic, size, timepoint_hr, dose) %>%
  summarise(
    Primary_Media_Valid = first(Primary_Media_Valid),
    PBS_Sensitivity_Valid = first(PBS_Sensitivity_Valid),
    Strict_Media_Valid = first(Strict_Media_Valid),
    Strict_PBS_Valid = first(Strict_PBS_Valid),
    Media_QC = first(Media_QC),
    PBS_QC = first(PBS_QC),
    Mean_Dead_Media = if_else(any(is.finite(Percent_Dead_Media)), mean(Percent_Dead_Media, na.rm = TRUE), NA_real_),
    SD_Tech_Media = if_else(sum(is.finite(Percent_Dead_Media)) >= 2, sd(Percent_Dead_Media, na.rm = TRUE), NA_real_),
    Mean_Dead_PBS = if_else(any(is.finite(Percent_Dead_PBS)), mean(Percent_Dead_PBS, na.rm = TRUE), NA_real_),
    SD_Tech_PBS = if_else(sum(is.finite(Percent_Dead_PBS)) >= 2, sd(Percent_Dead_PBS, na.rm = TRUE), NA_real_),
    Mean_Dead_MediaMedian = if_else(any(is.finite(Percent_Dead_MediaMedian)), mean(Percent_Dead_MediaMedian, na.rm = TRUE), NA_real_),
    SD_Tech_MediaMedian = if_else(sum(is.finite(Percent_Dead_MediaMedian)) >= 2, sd(Percent_Dead_MediaMedian, na.rm = TRUE), NA_real_),
    N_Tech_Media = sum(is.finite(Percent_Dead_Media)),
    N_Tech_PBS = sum(is.finite(Percent_Dead_PBS)),
    N_Tech_MediaMedian = sum(is.finite(Percent_Dead_MediaMedian)),
    .groups = "drop"
  ) %>%
  arrange(cell_line, replicate, microplastic, size, timepoint_hr, dose)

stopifnot(nrow(bio_dose_summary) == 1440)
write_csv(bio_dose_summary, file.path(out_dir, "MNP_v5_2_BiologicalReplicate_DoseMeans.csv"))

# =============================================================================
# 9. CC50 HELPERS: RETAIN ALL 480 DESIGN CELLS, FLAG QC FAILURES AS NA
# =============================================================================

estimate_group_cc50 <- function(df, response_col, valid_col = NULL, qc_fail_label = "QC_FAILED") {
  if (!is.null(valid_col)) {
    v <- unique(df[[valid_col]])
    v <- v[!is.na(v)]
    if (length(v) == 0 || !all(v)) {
      return(tibble(
        cc50 = NA_real_,
        cc50_censor = qc_fail_label,
        cc50_label = NA_character_,
        monotonic_nondecreasing = NA,
        crossing_interval = NA_character_
      ))
    }
  }
  estimate_cc50(df$dose, df[[response_col]])
}

# ---- Primary media --------------------------------------------------------------
cc50_media_all <- bio_dose_summary %>%
  group_by(cell_line, replicate, microplastic, size, timepoint_hr) %>%
  group_modify(~ estimate_group_cc50(.x, "Mean_Dead_Media", "Primary_Media_Valid", "QC_FAILED_MEDIA")) %>%
  ungroup() %>%
  mutate(background = "Media_primary_QC") %>%
  rename(Treatment = microplastic, Size = size, Timepoint = timepoint_hr, Replicate = replicate) %>%
  arrange(cell_line, Replicate, Treatment, Size, Timepoint)

stopifnot(nrow(cc50_media_all) == 480)
write_csv(cc50_media_all, file.path(out_dir, "CC50_ReplicateValues_PRIMARY_ALL480_QC.csv"))

cc50_primary_valid <- cc50_media_all %>%
  filter(is.finite(cc50), cc50_censor != "QC_FAILED_MEDIA")

# Compatibility table for downstream models: only QC-valid primary observations.
cc50_reps_compat <- cc50_primary_valid %>%
  transmute(
    cell_line,
    Size,
    Timepoint = as.character(Timepoint),
    Treatment,
    Replicate,
    cc50,
    cc50_censor,
    cc50_label,
    monotonic_nondecreasing,
    crossing_interval
  )

write_csv(cc50_reps_compat, file.path(out_dir, "CC50_ReplicateValues.csv"))
write_csv(cc50_reps_compat, file.path(out_dir, "CC50_ReplicateValues_PRIMARY_QCvalid.csv"))

# ---- Strict media-separation sensitivity (SSMD >= 1) ----------------------------
cc50_media_strict_all <- bio_dose_summary %>%
  group_by(cell_line, replicate, microplastic, size, timepoint_hr) %>%
  group_modify(~ estimate_group_cc50(.x, "Mean_Dead_Media", "Strict_Media_Valid", "QC_FAILED_STRICT_MEDIA")) %>%
  ungroup() %>%
  mutate(background = "Media_strict_SSMD_sensitivity") %>%
  rename(Treatment = microplastic, Size = size, Timepoint = timepoint_hr, Replicate = replicate) %>%
  arrange(cell_line, Replicate, Treatment, Size, Timepoint)

write_csv(cc50_media_strict_all, file.path(out_dir, "CC50_ReplicateValues_MEDIA_STRICT_SSMD_SENSITIVITY_ALL480.csv"))

# ---- Same-plate PBS sensitivity -------------------------------------------------
cc50_pbs_all <- bio_dose_summary %>%
  group_by(cell_line, replicate, microplastic, size, timepoint_hr) %>%
  group_modify(~ estimate_group_cc50(.x, "Mean_Dead_PBS", "PBS_Sensitivity_Valid", "QC_FAILED_PBS")) %>%
  ungroup() %>%
  mutate(background = "PBS_sensitivity") %>%
  rename(Treatment = microplastic, Size = size, Timepoint = timepoint_hr, Replicate = replicate) %>%
  arrange(cell_line, Replicate, Treatment, Size, Timepoint)

write_csv(cc50_pbs_all, file.path(out_dir, "CC50_ReplicateValues_PBS_sensitivity_ALL480.csv"))

# ---- Same-plate median-media sensitivity ---------------------------------------
# Median control values are intentionally sensitivity-only even if they would
# rescue the sign of the dynamic range.
cc50_media_median_all <- bio_dose_summary %>%
  group_by(cell_line, replicate, microplastic, size, timepoint_hr) %>%
  group_modify(~ estimate_cc50(.x$dose, .x$Mean_Dead_MediaMedian)) %>%
  ungroup() %>%
  mutate(background = "Media_median_sensitivity") %>%
  rename(Treatment = microplastic, Size = size, Timepoint = timepoint_hr, Replicate = replicate) %>%
  arrange(cell_line, Replicate, Treatment, Size, Timepoint)

write_csv(cc50_media_median_all, file.path(out_dir, "CC50_ReplicateValues_MediaMedian_sensitivity_ALL480.csv"))

# ---- Hierarchical same-plate rescue table (SENSITIVITY ONLY) --------------------
# Priority: primary media -> PBS -> median media. No controls are borrowed from a
# different plate. This table is not used as the manuscript primary analysis.
cc50_rescue <- cc50_media_all %>%
  select(cell_line, Replicate, Treatment, Size, Timepoint,
         media_cc50 = cc50, media_censor = cc50_censor) %>%
  left_join(
    cc50_pbs_all %>%
      select(cell_line, Replicate, Treatment, Size, Timepoint,
             pbs_cc50 = cc50, pbs_censor = cc50_censor),
    by = c("cell_line", "Replicate", "Treatment", "Size", "Timepoint")
  ) %>%
  left_join(
    cc50_media_median_all %>%
      select(cell_line, Replicate, Treatment, Size, Timepoint,
             median_media_cc50 = cc50, median_media_censor = cc50_censor),
    by = c("cell_line", "Replicate", "Treatment", "Size", "Timepoint")
  ) %>%
  mutate(
    Rescue_Source = case_when(
      is.finite(media_cc50) & media_censor != "QC_FAILED_MEDIA" ~ "Primary_Media",
      is.finite(pbs_cc50) & pbs_censor != "QC_FAILED_PBS" ~ "PBS_same_plate",
      is.finite(median_media_cc50) ~ "Median_Media_same_plate",
      TRUE ~ "No_same_plate_rescue"
    ),
    cc50_rescue = case_when(
      Rescue_Source == "Primary_Media" ~ media_cc50,
      Rescue_Source == "PBS_same_plate" ~ pbs_cc50,
      Rescue_Source == "Median_Media_same_plate" ~ median_media_cc50,
      TRUE ~ NA_real_
    )
  )

write_csv(cc50_rescue, file.path(out_dir, "CC50_ReplicateValues_SAMEPLATE_RESCUE_SENSITIVITY.csv"))

# =============================================================================
# 10. CONDITION-LEVEL PRIMARY CC50 SUMMARY
# =============================================================================

cc50_summary <- cc50_reps_compat %>%
  group_by(cell_line, Size, Timepoint, Treatment) %>%
  summarise(
    n = sum(is.finite(cc50)),
    n_expected = 3L,
    n_missing_qc = n_expected - n,
    cc50_mean = if_else(n > 0, mean(cc50, na.rm = TRUE), NA_real_),
    cc50_sd = if_else(n >= 2, sd(cc50, na.rm = TRUE), NA_real_),
    cc50_sem = if_else(n >= 2, cc50_sd / sqrt(n), NA_real_),
    n_lower_censored = sum(cc50_censor == "Lower", na.rm = TRUE),
    n_upper_censored = sum(cc50_censor == "Upper", na.rm = TRUE),
    n_unresolved = sum(cc50_censor == "Unresolved_nonmonotonic", na.rm = TRUE),
    n_nonmonotonic = sum(!monotonic_nondecreasing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cell_line, Treatment, Size, as.numeric(Timepoint))

# Keep all 160 condition rows even if a condition loses all three replicates.
all_conditions <- expand_grid(
  cell_line = CELL_LEVELS,
  Size = SIZE_LEVELS,
  Timepoint = as.character(TP_LEVELS),
  Treatment = POLY_LEVELS
)
cc50_summary <- all_conditions %>%
  left_join(cc50_summary, by = c("cell_line", "Size", "Timepoint", "Treatment")) %>%
  mutate(
    n = replace_na(n, 0L),
    n_expected = replace_na(n_expected, 3L),
    n_missing_qc = replace_na(n_missing_qc, 3L),
    n_lower_censored = replace_na(n_lower_censored, 0L),
    n_upper_censored = replace_na(n_upper_censored, 0L),
    n_unresolved = replace_na(n_unresolved, 0L),
    n_nonmonotonic = replace_na(n_nonmonotonic, 0L)
  ) %>%
  arrange(match(cell_line, CELL_LEVELS), match(Treatment, POLY_LEVELS),
          match(Size, SIZE_LEVELS), as.numeric(Timepoint))

stopifnot(nrow(cc50_summary) == 160)
write_csv(cc50_summary, file.path(out_dir, "Cell_Death_CC50_Summary_withErrorBars_REGEN.csv"))

# =============================================================================
# 11. QC SENSITIVITY SUMMARIES
# =============================================================================

plate_qc_counts <- control_summary %>%
  summarise(
    total_plates = n(),
    primary_media_valid = sum(Primary_Media_Valid),
    primary_media_failed = sum(!Primary_Media_Valid),
    media_pass = sum(Media_QC == "PASS"),
    media_warn = sum(Media_QC == "WARN_separation"),
    media_fail_direction = sum(Media_QC == "FAIL_direction"),
    media_fail_separation = sum(Media_QC == "FAIL_separation"),
    pbs_same_plate_rescue_candidates = sum(PBS_Rescue_Candidate),
    median_media_rescue_candidates = sum(Median_Media_Rescue_Candidate)
  )
write_csv(plate_qc_counts, file.path(out_dir, "MNP_v5_2_PlateQC_OverallSummary.csv"))

cc50_qc_summary <- cc50_media_all %>%
  group_by(cell_line) %>%
  summarise(
    n_design_cc50 = n(),
    n_qc_failed = sum(cc50_censor == "QC_FAILED_MEDIA", na.rm = TRUE),
    n_valid_cc50 = sum(is.finite(cc50)),
    n_interpolated = sum(cc50_censor == "None", na.rm = TRUE),
    n_lower = sum(cc50_censor == "Lower", na.rm = TRUE),
    n_upper = sum(cc50_censor == "Upper", na.rm = TRUE),
    n_unresolved = sum(cc50_censor == "Unresolved_nonmonotonic", na.rm = TRUE),
    n_nonmonotonic_valid = sum(!monotonic_nondecreasing & is.finite(cc50), na.rm = TRUE),
    .groups = "drop"
  )
write_csv(cc50_qc_summary, file.path(out_dir, "MNP_v5_2_CC50_QC_byCellLine.csv"))

rescue_summary <- cc50_rescue %>%
  count(cell_line, Rescue_Source, name = "n_cc50") %>%
  arrange(cell_line, Rescue_Source)
write_csv(rescue_summary, file.path(out_dir, "MNP_v5_2_SamePlateRescue_Summary.csv"))

# =============================================================================
# 12. QC PLOTS
# =============================================================================

# Plate control separation plot
p_controls <- control_summary %>%
  mutate(
    cell_line = factor(cell_line, levels = CELL_LEVELS),
    replicate = factor(replicate),
    timepoint_hr = factor(timepoint_hr, levels = TP_LEVELS)
  ) %>%
  ggplot(aes(x = timepoint_hr, y = SSMD_Media, shape = replicate)) +
  geom_hline(yintercept = STRICT_SSMD_MIN, linetype = 2) +
  geom_hline(yintercept = WARN_SSMD_MIN, linetype = 3) +
  geom_point(size = 2.2) +
  facet_wrap(~ cell_line, scales = "free_y") +
  labs(
    x = "Timepoint (h)", y = "Signed SSMD: Methanol vs media",
    shape = "Biological replicate",
    title = "Plate-control separation QC",
    subtitle = paste0("Primary exclusion only when Methanol <= Media; SSMD < ", STRICT_SSMD_MIN,
                      " is strict-sensitivity flag; ", STRICT_SSMD_MIN, " to <", WARN_SSMD_MIN, " is weak-separation flag")
  ) +
  theme_bw(base_size = 9)

ggsave(file.path(out_dir, "MNP_v5_2_PlateControl_SSMD_QC.png"),
       p_controls, width = 10, height = 6, dpi = 300)

# Dose-response plot: only QC-valid primary media plates
plot_df <- bio_dose_summary %>%
  filter(Primary_Media_Valid, is.finite(Mean_Dead_Media)) %>%
  mutate(
    cell_line = factor(cell_line, levels = CELL_LEVELS),
    microplastic = factor(microplastic, levels = POLY_LEVELS),
    size = factor(size, levels = SIZE_LEVELS),
    timepoint_hr = factor(timepoint_hr, levels = TP_LEVELS),
    replicate = factor(replicate)
  )

p_dose <- ggplot(
  plot_df,
  aes(x = dose, y = Mean_Dead_Media,
      group = interaction(microplastic, replicate), linetype = replicate)
) +
  geom_hline(yintercept = 50, linetype = 3) +
  geom_line(aes(color = microplastic), linewidth = 0.45, alpha = 0.8) +
  geom_point(aes(color = microplastic), size = 1.2) +
  scale_x_log10(breaks = DOSE_LEVELS, labels = DOSE_LEVELS) +
  facet_grid(size + cell_line ~ timepoint_hr) +
  labs(
    x = "MNP concentration (ug/mL; log scale)",
    y = "Mean cell death within biological replicate (%)",
    color = "Polymer", linetype = "Biological replicate",
    title = "MNP dose-response QC: QC-valid media-normalized plates"
  ) +
  theme_bw(base_size = 9)

ggsave(file.path(out_dir, "MNP_v5_2_DoseResponse_PRIMARY_QCvalid.png"),
       p_dose, width = 14, height = 15, dpi = 300)

# =============================================================================
# 13. FINAL AUDIT
# =============================================================================

cat("\n==================== v5.2 COMPLETE ====================\n")
cat("Input directory: ", raw_dir, "\n", sep = "")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Raw MARS plates parsed: ", length(mars_files), " / 60 expected\n", sep = "")
cat("Corrected doses: ", paste(DOSE_LEVELS, collapse = ", "), " ug/mL\n", sep = "")
cat("Raw experimental wells: ", nrow(samples_raw), " / 4320 expected\n", sep = "")
cat("Raw control wells: ", nrow(controls_raw), " / 540 expected\n", sep = "")
cat("Biological-replicate dose means: ", nrow(bio_dose_summary), " / 1440 expected\n", sep = "")
cat("Design-level replicate CC50 cells: ", nrow(cc50_media_all), " / 480 expected\n", sep = "")
cat("QC-valid primary media CC50s: ", nrow(cc50_reps_compat), "\n", sep = "")
cat("Strict SSMD>=1 sensitivity CC50s: ", sum(is.finite(cc50_media_strict_all$cc50)), "\n", sep = "")
cat("Primary media plates valid: ", sum(control_summary$Primary_Media_Valid), " / 60\n", sep = "")
cat("Primary media plates failed: ", sum(!control_summary$Primary_Media_Valid), " / 60\n", sep = "")
cat("  PASS: ", sum(control_summary$Media_QC == "PASS"), "\n", sep = "")
cat("  WARN separation (1 to <2): ", sum(control_summary$Media_QC == "WARN_separation"), "\n", sep = "")
cat("  FAIL direction: ", sum(control_summary$Media_QC == "FAIL_direction"), "\n", sep = "")
cat("  WARN poor separation (<1): ", sum(control_summary$Media_QC == "WARN_poor_separation"), "\n", sep = "")
cat("PBS same-plate rescue candidate plates: ", sum(control_summary$PBS_Rescue_Candidate), "\n", sep = "")
cat("Median-media same-plate rescue candidate plates: ", sum(control_summary$Median_Media_Rescue_Candidate), "\n", sep = "")
cat("\nPrimary CC50 QC by cell line:\n")
print(cc50_qc_summary)
cat("\nSame-plate rescue hierarchy (SENSITIVITY ONLY):\n")
print(rescue_summary, n = Inf)
cat("\nTechnical-replicate QC flags:\n")
print(tech_qc %>% count(Technical_QC_Flag, name = "n_groups") %>% arrange(desc(n_groups)), n = Inf)
cat("=======================================================\n")

if (any(!control_summary$Primary_Media_Valid)) {
  warning(
    "Some plates failed primary media-control QC. Their primary media-normalized values and CC50s were set to NA. ",
    "Review MNP_v5_2_PlateControl_QC.csv and sensitivity outputs; no cross-plate control borrowing was performed."
  )
}
