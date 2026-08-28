# ==============================================================================
# 09_make_supplementary_table_s2.R
#
# PURPOSE
#   Build Supplementary Table S2 describing the no-cell particle-fluorescence
#   control and the final interference diagnostics.
#
# OUTPUT
#   manuscript_outputs/Supplementary_Table_S2_ParticleInterference.xlsx
#
# SHEETS
#   S2A_NoCellSignal
#   S2B_Correlations
#   S2C_Regression
#
# IMPORTANT
#   This workbook intentionally does NOT include the repeated-observation mixed
#   model from v6.5 because the particle-only predictor was measured at only 24
#   unique polymer x size x dose conditions. Final inferential emphasis is
#   placed on condition-level analyses with n = 24.
# ==============================================================================

# Run this script from the repository root.
project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(project_dir, "scripts")) ||
    !dir.exists(file.path(project_dir, "results"))) {
  stop("Run this script from the repository root.")
}


suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

v65_dir <- file.path(project_dir, "results", "interference", "v6_5")
v66_dir <- file.path(project_dir, "results", "interference", "v6_6")
out_dir <- file.path(project_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(
  out_dir,
  "Supplementary_Table_S2_ParticleInterference.xlsx"
)

read_req <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

# ------------------------------------------------------------------------------
# S2A. NO-CELL PARTICLE-ASSOCIATED FLUORESCENCE
# ------------------------------------------------------------------------------

s2a <- read_req(
  file.path(v65_dir, "V6_5_NoCell_ParticleFluorescence_ByCondition.csv")
) %>%
  transmute(
    Polymer = microplastic,
    `Particle-size class` = size,
    `MNP concentration (µg/mL)` = dose,
    `No-cell replicates` = n_no_cell,
    `Mean no-cell fluorescence (FI)` = no_cell_mean,
    `SD no-cell fluorescence (FI)` = no_cell_sd,
    `No-cell media mean (FI)` = no_cell_media_mean,
    `Particle-associated excess fluorescence (FI)` = particle_excess_signed,
    `Fold versus no-cell media` = fold_vs_no_cell_media
  ) %>%
  arrange(
    factor(Polymer, levels = c("PE", "PS", "PVC", "Mix")),
    factor(`Particle-size class`, levels = c("Small", "Large")),
    `MNP concentration (µg/mL)`
  )

# ------------------------------------------------------------------------------
# S2B. CORRELATION DIAGNOSTICS
# ------------------------------------------------------------------------------

overall <- read_req(
  file.path(v65_dir, "V6_5_Spearman_Overall.csv")
) %>%
  transmute(
    Analysis = "Overall",
    Stratum = outcome,
    `Independent conditions (n)` = n_conditions,
    `Spearman rho` = spearman_rho,
    `P value` = p_value,
    `BH q value` = q_BH
  )

by_cell <- read_req(
  file.path(v65_dir, "V6_5_Spearman_ByCellLine.csv")
) %>%
  mutate(
    cell_line = recode(
      cell_line,
      "HTR8" = "HTR-8/SVneo",
      "JEG3" = "JEG-3",
      .default = cell_line
    )
  ) %>%
  transmute(
    Analysis = "By cell line",
    Stratum = cell_line,
    `Independent conditions (n)` = n_conditions,
    `Spearman rho` = spearman_rho,
    `P value` = p_value,
    `BH q value` = q_BH
  )

by_time <- read_req(
  file.path(v65_dir, "V6_5_Spearman_ByTimepoint.csv")
) %>%
  transmute(
    Analysis = "By exposure duration",
    Stratum = paste0(timepoint_hr, " h"),
    `Independent conditions (n)` = n_conditions,
    `Spearman rho` = spearman_rho,
    `P value` = p_value,
    `BH q value` = q_BH
  )

within_dose <- read_req(
  file.path(v65_dir, "V6_5_Spearman_WithinDose.csv")
) %>%
  transmute(
    Analysis = "Within MNP concentration",
    Stratum = paste0(dose, " µg/mL"),
    `Independent conditions (n)` = n_conditions,
    `Spearman rho` = spearman_rho,
    `P value` = p_value,
    `BH q value` = q_BH
  )

s2b <- bind_rows(overall, by_cell, by_time, within_dose)

# ------------------------------------------------------------------------------
# S2C. FINAL CONDITION-LEVEL REGRESSION
# ------------------------------------------------------------------------------

effects <- read_req(
  file.path(v66_dir, "V6_6_ParticleSignal_EffectSummary.csv")
) %>%
  transmute(
    Analysis = "Condition-level linear regression",
    Outcome = outcome,
    `Independent conditions (n)` = n_conditions,
    `Effect metric` = "Change in apparent cell death (%) per 1-SD increase in particle-only fluorescence",
    Estimate = estimate_percent_death_per_1SD_particle_excess,
    `95% CI lower` = ci_lower,
    `95% CI upper` = ci_upper,
    `P value` = p_value,
    `BH q value` = q_BH
  )

hc3_file <- file.path(v66_dir, "V6_6_ParticleSignal_HC3_Sensitivity.csv")

if (file.exists(hc3_file)) {
  hc3 <- read_req(hc3_file) %>%
    transmute(
      Analysis = "HC3 robust-SE sensitivity",
      Outcome = outcome,
      `Independent conditions (n)` = 24L,
      `Effect metric` = "Change in apparent cell death (%) per 1-SD increase in particle-only fluorescence",
      Estimate = estimate,
      `95% CI lower` = NA_real_,
      `95% CI upper` = NA_real_,
      `P value` = p_value,
      `BH q value` = q_BH
    )
  s2c <- bind_rows(effects, hc3)
} else {
  s2c <- effects
}

# ------------------------------------------------------------------------------
# WORKBOOK
# ------------------------------------------------------------------------------

wb <- createWorkbook()

title_style <- createStyle(
  fgFill = "#1F4E78",
  fontColour = "#FFFFFF",
  textDecoration = "bold",
  fontSize = 12,
  valign = "center"
)

note_style <- createStyle(
  fgFill = "#F3F6F8",
  textDecoration = "italic",
  fontColour = "#404040",
  fontSize = 9,
  wrapText = TRUE,
  valign = "center"
)

header_style <- createStyle(
  fgFill = "#D9EAF7",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "Bottom",
  borderColour = "#7F8C8D"
)

sig_style <- createStyle(
  fgFill = "#E2F0D9",
  textDecoration = "bold"
)

add_table_sheet <- function(
  wb, sheet_name, title, note, dat,
  col_widths = NULL
) {
  addWorksheet(wb, sheet_name)

  writeData(wb, sheet_name, title, startRow = 1, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:ncol(dat), rows = 1)
  addStyle(
    wb, sheet_name, title_style,
    rows = 1, cols = 1:ncol(dat), gridExpand = TRUE
  )
  setRowHeights(wb, sheet_name, rows = 1, heights = 26)

  writeData(wb, sheet_name, note, startRow = 2, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:ncol(dat), rows = 2)
  addStyle(
    wb, sheet_name, note_style,
    rows = 2, cols = 1:ncol(dat), gridExpand = TRUE
  )
  setRowHeights(wb, sheet_name, rows = 2, heights = 42)

  writeData(
    wb, sheet_name, dat,
    startRow = 4, startCol = 1,
    headerStyle = header_style
  )

  freezePane(wb, sheet_name, firstActiveRow = 5)

  if (is.null(col_widths)) {
    setColWidths(wb, sheet_name, cols = 1:ncol(dat), widths = "auto")
  } else {
    setColWidths(
      wb, sheet_name,
      cols = seq_along(col_widths),
      widths = col_widths
    )
  }

  # P/q formatting
  pcols <- which(names(dat) %in% c("P value", "BH q value"))
  if (length(pcols) > 0) {
    addStyle(
      wb, sheet_name,
      createStyle(numFmt = "0.0000"),
      rows = 5:(4 + nrow(dat)),
      cols = pcols,
      gridExpand = TRUE
    )
  }

  qcol <- which(names(dat) == "BH q value")
  if (length(qcol) == 1) {
    conditionalFormatting(
      wb, sheet_name,
      cols = qcol,
      rows = 5:(4 + nrow(dat)),
      rule = "<0.05",
      style = sig_style
    )
  }
}

add_table_sheet(
  wb,
  "S2A_NoCellSignal",
  "Supplementary Table S2A. Particle-associated fluorescence in no-cell controls",
  paste0(
    "No-cell fluorescence was measured for each polymer × particle-size × concentration condition ",
    "using the EthD-1 acquisition settings. Particle-associated excess fluorescence was calculated ",
    "as the condition-specific no-cell MNP signal minus the mean no-cell media signal. These values ",
    "were used to characterize potential optical interference and were not directly subtracted from ",
    "cell-containing wells."
  ),
  s2a,
  col_widths = c(12, 18, 20, 17, 24, 22, 21, 30, 22)
)

add_table_sheet(
  wb,
  "S2B_Correlations",
  "Supplementary Table S2B. Association between particle-only fluorescence and apparent EthD-1-defined cell death",
  paste0(
    "Spearman correlations evaluate whether particle-associated fluorescence tracked the original ",
    "apparent EthD-1-defined cell-death signal across unique polymer × particle-size × concentration ",
    "conditions. Analyses were performed overall and stratified by cell line, exposure duration, and MNP concentration."
  ),
  s2b,
  col_widths = c(24, 28, 22, 18, 16, 16)
)

add_table_sheet(
  wb,
  "S2C_Regression",
  "Supplementary Table S2C. Condition-level regression of particle-only fluorescence and apparent cell death",
  paste0(
    "Linear regression models used the 24 unique polymer × particle-size × concentration conditions ",
    "and adjusted for polymer identity, particle-size class, and concentration. Estimates represent the ",
    "change in apparent EthD-1-defined cell death per 1-SD increase in particle-only fluorescence. ",
    "HC3 robust-standard-error sensitivity results are shown when available."
  ),
  s2c,
  col_widths = c(26, 28, 22, 48, 15, 15, 15, 16, 16)
)

saveWorkbook(wb, out_file, overwrite = TRUE)

cat("\n================ SUPPLEMENTARY TABLE S2 COMPLETE ================\n")
cat("Output: ", out_file, "\n", sep = "")
cat("S2A rows: ", nrow(s2a), "\n", sep = "")
cat("S2B rows: ", nrow(s2b), "\n", sep = "")
cat("S2C rows: ", nrow(s2c), "\n", sep = "")
cat("=================================================================\n")
