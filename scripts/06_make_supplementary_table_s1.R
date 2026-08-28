# ==============================================================================
# 06_make_supplementary_table_s1.R
# Generate Supplementary Table S1 workbook from final v6.3 CC50 outputs.
#
# Output:
#   manuscript_outputs/Supplementary_Table_S1_CC50_v6_3.xlsx
#
# Sheets:
#   S1A_Global
#   S1B_Cell_vs_BeWo
#   S1C_Exposure
#   S1D_Interactions
#   S1E_Int_Contrasts
#   S1F_AllPairwise
#   S1G_CellSummary
#   S1H_Condition
# ==============================================================================

# Run this script from the repository root.
project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(project_dir, "scripts")) ||
    !dir.exists(file.path(project_dir, "metadata"))) {
  stop("Run this script from the repository root (the directory containing scripts/ and metadata/).")
}


suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

res_dir <- file.path(project_dir, "results", "v6_3")
out_dir <- file.path(project_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(
  out_dir,
  "Supplementary_Table_S1_CC50_v6_3.xlsx"
)

display_cell <- function(x) {
  x <- str_replace_all(x, fixed("(THP-1)"), "THP-1")
  dplyr::recode(
    x,
    "HTR8" = "HTR-8/SVneo",
    "JEG3" = "JEG-3",
    .default = x
  )
}

display_contrast <- function(x) {
  x %>%
    str_replace_all(fixed("(THP-1)"), "THP-1") %>%
    str_replace_all("HTR8", "HTR-8/SVneo") %>%
    str_replace_all("JEG3", "JEG-3")
}

display_term <- function(x) {
  recode(
    x,
    "cell_line" = "Cell line",
    "Treatment" = "Polymer treatment",
    "Size" = "Particle-size class",
    "Timepoint" = "Exposure duration",
    "cell_line:Treatment" = "Cell line × polymer treatment",
    "cell_line:Size" = "Cell line × particle-size class",
    "cell_line:Timepoint" = "Cell line × exposure duration",
    .default = x
  )
}

read_res <- function(filename) {
  fp <- file.path(res_dir, filename)
  if (!file.exists(fp)) stop("Missing v6.3 output: ", fp)
  readr::read_csv(fp, show_col_types = FALSE)
}

# ------------------------------------------------------------------------------
# S1A: Global tests
# ------------------------------------------------------------------------------

s1a_primary <- read_res("V6_3_PRIMARY_Additive_LMM_TypeIII.csv") %>%
  transmute(
    Model = "Primary additive LMM",
    Term = display_term(Term),
    Statistic = "F",
    `Statistic value` = F_value,
    `Numerator df` = NumDF,
    `Denominator df` = DenDF,
    `P value` = p,
    `BH q value` = q_BH
  )

s1a_ordinal <- read_res("V6_3_SENS_OrdinalCLMM_LRT_Tests.csv") %>%
  transmute(
    Model = "Censoring-aware ordinal CLMM",
    Term = display_term(Term),
    Statistic = "Likelihood-ratio χ²",
    `Statistic value` = LRT,
    `Numerator df` = NA_real_,
    `Denominator df` = NA_real_,
    `P value` = p,
    `BH q value` = q_BH
  )

s1a_strict <- read_res("V6_3_SENS_StrictSSMD_Additive_TypeIII.csv") %>%
  transmute(
    Model = "Strict SSMD ≥1 additive LMM",
    Term = display_term(Term),
    Statistic = "F",
    `Statistic value` = F_value,
    `Numerator df` = NumDF,
    `Denominator df` = DenDF,
    `P value` = p,
    `BH q value` = q_BH
  )

s1a <- bind_rows(s1a_primary, s1a_ordinal, s1a_strict)

# ------------------------------------------------------------------------------
# S1B: Cell line versus BeWo
# ------------------------------------------------------------------------------

s1b_primary <- read_res("V6_3_PRIMARY_CellLine_vs_BeWo.csv") %>%
  transmute(
    Model = "Primary additive LMM",
    Contrast = paste0(display_cell(test_cell), " vs BeWo"),
    `Effect metric` = "ΔCC50 (µg/mL)",
    Estimate,
    `95% CI lower` = lowerCL,
    `95% CI upper` = upperCL,
    `P value` = p_raw,
    `BH q value` = q_BH
  )

s1b_ordinal <- read_res("V6_3_SENS_OrdinalCLMM_CellLine_vs_BeWo.csv") %>%
  transmute(
    Model = "Censoring-aware ordinal CLMM",
    Contrast = paste0(display_cell(test_cell), " vs BeWo"),
    `Effect metric` = "OR for higher CC50 category",
    Estimate = OR_higher_CC50_category,
    `95% CI lower` = OR_lower95,
    `95% CI upper` = OR_upper95,
    `P value` = p_raw,
    `BH q value` = q_BH
  )

s1b_strict <- read_res("V6_3_SENS_StrictSSMD_CellLine_vs_BeWo.csv") %>%
  transmute(
    Model = "Strict SSMD ≥1 additive LMM",
    Contrast = paste0(display_cell(test_cell), " vs BeWo"),
    `Effect metric` = "ΔCC50 (µg/mL)",
    Estimate,
    `95% CI lower` = lowerCL,
    `95% CI upper` = upperCL,
    `P value` = p_raw,
    `BH q value` = q_BH
  )

s1b <- bind_rows(s1b_primary, s1b_ordinal, s1b_strict)

# ------------------------------------------------------------------------------
# S1C: Exposure-factor contrasts
# ------------------------------------------------------------------------------

primary_exposure <- bind_rows(
  read_res("V6_3_PRIMARY_Treatment_Pairwise.csv") %>% mutate(Factor = "Polymer treatment"),
  read_res("V6_3_PRIMARY_Size_Pairwise.csv") %>% mutate(Factor = "Particle-size class"),
  read_res("V6_3_PRIMARY_Timepoint_Pairwise.csv") %>% mutate(Factor = "Exposure duration")
) %>%
  transmute(
    Model = "Primary additive LMM",
    Factor,
    Contrast = contrast,
    `Effect metric` = "ΔCC50 (µg/mL)",
    Estimate,
    `95% CI lower` = lowerCL,
    `95% CI upper` = upperCL,
    `P value` = p_raw,
    `BH q value` = q_BH
  )

ordinal_exposure <- bind_rows(
  read_res("V6_3_SENS_OrdinalCLMM_Treatment_Pairwise.csv") %>% mutate(Factor = "Polymer treatment"),
  read_res("V6_3_SENS_OrdinalCLMM_Size_Pairwise.csv") %>% mutate(Factor = "Particle-size class"),
  read_res("V6_3_SENS_OrdinalCLMM_Timepoint_Pairwise.csv") %>% mutate(Factor = "Exposure duration")
) %>%
  transmute(
    Model = "Censoring-aware ordinal CLMM",
    Factor,
    Contrast = contrast,
    `Effect metric` = "OR for higher CC50 category",
    Estimate = OR_higher_CC50_category,
    `95% CI lower` = OR_lower95,
    `95% CI upper` = OR_upper95,
    `P value` = p_raw,
    `BH q value` = q_BH
  )

s1c <- bind_rows(primary_exposure, ordinal_exposure)

# ------------------------------------------------------------------------------
# S1D: Secondary interaction omnibus tests
# ------------------------------------------------------------------------------

s1d <- bind_rows(
  read_res("V6_3_SECONDARY_CellLine_x_Treatment_TypeIII.csv") %>%
    mutate(`Secondary model` = "Cell line × polymer treatment"),
  read_res("V6_3_SECONDARY_CellLine_x_Timepoint_TypeIII.csv") %>%
    mutate(`Secondary model` = "Cell line × exposure duration"),
  read_res("V6_3_SECONDARY_CellLine_x_Size_TypeIII.csv") %>%
    mutate(`Secondary model` = "Cell line × particle-size class")
) %>%
  transmute(
    `Secondary model`,
    Term = display_term(Term),
    `Numerator df` = NumDF,
    `Denominator df` = DenDF,
    `F value` = F_value,
    `P value` = p,
    `BH q value` = q_BH
  )

# ------------------------------------------------------------------------------
# S1E: Secondary interaction-specific BeWo contrasts
# ------------------------------------------------------------------------------

s1e_treatment <- read_res("V6_3_SECONDARY_CellLine_x_Treatment_BeWoContrasts.csv") %>%
  transmute(
    `Secondary model` = "Cell line × polymer treatment",
    Modifier = "Polymer treatment",
    Level = Treatment,
    Contrast = paste0(display_cell(test_cell), " vs BeWo"),
    `ΔCC50 (µg/mL)` = Estimate,
    `95% CI lower` = lowerCL,
    `95% CI upper` = upperCL,
    `P value` = p_raw,
    `BH q within panel` = q_BH_within_panel
  )

s1e_time <- read_res("V6_3_SECONDARY_CellLine_x_Timepoint_BeWoContrasts.csv") %>%
  transmute(
    `Secondary model` = "Cell line × exposure duration",
    Modifier = "Exposure duration",
    Level = paste0(Timepoint, " h"),
    Contrast = paste0(display_cell(test_cell), " vs BeWo"),
    `ΔCC50 (µg/mL)` = Estimate,
    `95% CI lower` = lowerCL,
    `95% CI upper` = upperCL,
    `P value` = p_raw,
    `BH q within panel` = q_BH_within_panel
  )

s1e_size <- read_res("V6_3_SECONDARY_CellLine_x_Size_BeWoContrasts.csv") %>%
  transmute(
    `Secondary model` = "Cell line × particle-size class",
    Modifier = "Particle-size class",
    Level = Size,
    Contrast = paste0(display_cell(test_cell), " vs BeWo"),
    `ΔCC50 (µg/mL)` = Estimate,
    `95% CI lower` = lowerCL,
    `95% CI upper` = upperCL,
    `P value` = p_raw,
    `BH q within panel` = q_BH_within_panel
  )

s1e <- bind_rows(s1e_treatment, s1e_time, s1e_size)

# ------------------------------------------------------------------------------
# S1F: All pairwise cell-line comparisons
# ------------------------------------------------------------------------------

s1f_primary <- read_res("V6_3_SUPP_AllPairwise_CellLine_Contrasts.csv") %>%
  transmute(
    Model = "Primary additive LMM",
    Contrast = display_contrast(contrast),
    `Effect metric` = "ΔCC50 (µg/mL)",
    Estimate,
    `95% CI lower` = lowerCL,
    `95% CI upper` = upperCL,
    `P value` = p_raw,
    `BH q value` = q_BH
  )

s1f_ordinal <- read_res("V6_3_SENS_OrdinalCLMM_AllPairwise_CellLine.csv") %>%
  transmute(
    Model = "Censoring-aware ordinal CLMM",
    Contrast = display_contrast(contrast),
    `Effect metric` = "OR for higher CC50 category",
    Estimate = OR_higher_CC50_category,
    `95% CI lower` = OR_lower95,
    `95% CI upper` = OR_upper95,
    `P value` = p_raw,
    `BH q value` = q_BH
  )

s1f_strict <- read_res("V6_3_SENS_StrictSSMD_AllPairwise_CellLine.csv") %>%
  transmute(
    Model = "Strict SSMD ≥1 additive LMM",
    Contrast = display_contrast(contrast),
    `Effect metric` = "ΔCC50 (µg/mL)",
    Estimate,
    `95% CI lower` = lowerCL,
    `95% CI upper` = upperCL,
    `P value` = p_raw,
    `BH q value` = q_BH
  )

s1f <- bind_rows(s1f_primary, s1f_ordinal, s1f_strict)

# ------------------------------------------------------------------------------
# S1G/S1H: Descriptive summaries
# ------------------------------------------------------------------------------

s1g <- read_res("V6_3_CC50_ByCellLine_Descriptive.csv") %>%
  transmute(
    `Cell line` = display_cell(cell_line),
    `N CC50 estimates` = n,
    `N biological replicates` = n_bioreps,
    Mean = mean_boundary,
    SD = sd_boundary,
    Median = median_boundary,
    Q1 = q25_boundary,
    Q3 = q75_boundary,
    `N ≤0.2` = n_lower,
    `N interpolated` = n_interpolated,
    `N >20` = n_upper,
    `Proportion ≤0.2` = prop_lower,
    `Proportion interpolated` = prop_interpolated,
    `Proportion >20` = prop_upper
  )

s1h <- read_res("V6_3_CC50_ByCondition_Descriptive.csv") %>%
  transmute(
    `Cell line` = display_cell(cell_line),
    Polymer = Treatment,
    `Particle size` = Size,
    `Exposure duration` = paste0(Timepoint, " h"),
    N = n,
    `N biological replicates` = n_bioreps,
    Mean = mean_boundary,
    SD = sd_boundary,
    Median = median_boundary,
    Q1 = q25_boundary,
    Q3 = q75_boundary,
    `N ≤0.2` = n_lower,
    `N interpolated` = n_interpolated,
    `N >20` = n_upper,
    `Proportion ≤0.2` = prop_lower,
    `Proportion interpolated` = prop_interpolated,
    `Proportion >20` = prop_upper
  )

# ------------------------------------------------------------------------------
# Workbook
# ------------------------------------------------------------------------------

tables <- list(
  S1A_Global = list(
    title = "Supplementary Table S1A. Global tests from the primary and sensitivity models",
    note = "Primary additive LMM, censoring-aware ordinal CLMM, and strict SSMD≥1 sensitivity analysis. BH-adjusted q values are reported within each model.",
    data = s1a
  ),
  S1B_Cell_vs_BeWo = list(
    title = "Supplementary Table S1B. Adjusted cell-line comparisons relative to BeWo",
    note = "Positive ΔCC50 indicates higher boundary-assigned CC50 than BeWo. For the ordinal model, OR>1 indicates greater odds of occupying a higher CC50 category (≤0.2, 0.2–20, >20 µg/mL).",
    data = s1b
  ),
  S1C_Exposure = list(
    title = "Supplementary Table S1C. Adjusted polymer, particle-size, and exposure-duration contrasts",
    note = "Pairwise contrasts from the primary additive LMM and censoring-aware ordinal sensitivity model. Positive ΔCC50 or OR>1 indicates higher CC50 / lower cytotoxic susceptibility for the first level in the contrast.",
    data = s1c
  ),
  S1D_Interactions = list(
    title = "Supplementary Table S1D. Secondary interaction-model omnibus tests",
    note = "Secondary models evaluated one cell-line interaction at a time. These analyses are secondary; no cell-line interaction was significant after FDR correction.",
    data = s1d
  ),
  S1E_Int_Contrasts = list(
    title = "Supplementary Table S1E. Secondary interaction-specific cell-line contrasts relative to BeWo",
    note = "BeWo-referenced contrasts from the secondary interaction models. q values are BH-adjusted within the corresponding interaction panel and should be interpreted as secondary/exploratory.",
    data = s1e
  ),
  S1F_AllPairwise = list(
    title = "Supplementary Table S1F. All pairwise cell-line comparisons",
    note = "All pairwise cell-line contrasts from the primary additive LMM, censoring-aware ordinal CLMM, and strict SSMD≥1 sensitivity analysis.",
    data = s1f
  ),
  S1G_CellSummary = list(
    title = "Supplementary Table S1G. Descriptive boundary-assigned CC50 summary by cell line",
    note = "Descriptive summaries across all QC-valid conditions. Boundary values of 0.2 and 20 µg/mL represent lower- and upper-censored estimates, respectively.",
    data = s1g
  ),
  S1H_Condition = list(
    title = "Supplementary Table S1H. Condition-level descriptive and censoring summary",
    note = "Condition-level summaries used to construct Table 1 and Figure 2. Boundary-assigned values are descriptive; explicit censor counts and proportions are provided.",
    data = s1h
  )
)

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

# Contents
addWorksheet(wb, "Contents")
writeData(
  wb, "Contents",
  "Supplementary Table S1. Final CC50 statistical and descriptive results",
  startRow = 1, startCol = 1
)
mergeCells(wb, "Contents", cols = 1:3, rows = 1)
addStyle(wb, "Contents", title_style, rows = 1, cols = 1:3, gridExpand = TRUE)

contents <- tibble(
  Sheet = names(tables),
  Subtable = str_extract(names(tables), "^S1[A-H]"),
  Contents = map_chr(tables, ~str_remove(.x$title, "^Supplementary Table S1[A-H]\\. "))
)

writeData(wb, "Contents", contents, startRow = 3, headerStyle = header_style)
setColWidths(wb, "Contents", cols = 1, widths = 22)
setColWidths(wb, "Contents", cols = 2, widths = 10)
setColWidths(wb, "Contents", cols = 3, widths = 58)
freezePane(wb, "Contents", firstActiveRow = 4)

for (nm in names(tables)) {
  tab <- tables[[nm]]
  dat <- tab$data
  addWorksheet(wb, nm)

  writeData(wb, nm, tab$title, startRow = 1, startCol = 1)
  mergeCells(wb, nm, cols = 1:ncol(dat), rows = 1)
  addStyle(wb, nm, title_style, rows = 1, cols = 1:ncol(dat), gridExpand = TRUE)
  setRowHeights(wb, nm, rows = 1, heights = 24)

  writeData(wb, nm, tab$note, startRow = 2, startCol = 1)
  mergeCells(wb, nm, cols = 1:ncol(dat), rows = 2)
  addStyle(wb, nm, note_style, rows = 2, cols = 1:ncol(dat), gridExpand = TRUE)
  setRowHeights(wb, nm, rows = 2, heights = 38)

  writeData(wb, nm, dat, startRow = 4, headerStyle = header_style)
  freezePane(wb, nm, firstActiveRow = 5)

  setColWidths(wb, nm, cols = 1:ncol(dat), widths = 14)
  setColWidths(wb, nm, cols = 1, widths = 24)
  if (ncol(dat) >= 2) setColWidths(wb, nm, cols = 2, widths = 24)
  if (ncol(dat) >= 3) setColWidths(wb, nm, cols = 3, widths = 22)

  q_cols <- which(names(dat) %in% c("BH q value", "BH q within panel"))
  if (length(q_cols) > 0) {
    for (qc in q_cols) {
      conditionalFormatting(
        wb, nm,
        cols = qc,
        rows = 5:(4 + nrow(dat)),
        rule = "<0.05",
        style = sig_style
      )
    }
  }

  prop_cols <- which(str_detect(names(dat), "^Proportion "))
  if (length(prop_cols) > 0) {
    addStyle(
      wb, nm,
      createStyle(numFmt = "0.0%"),
      rows = 5:(4 + nrow(dat)),
      cols = prop_cols,
      gridExpand = TRUE
    )
  }

  p_cols <- which(names(dat) %in% c("P value", "BH q value", "BH q within panel"))
  if (length(p_cols) > 0) {
    addStyle(
      wb, nm,
      createStyle(numFmt = "0.0000"),
      rows = 5:(4 + nrow(dat)),
      cols = p_cols,
      gridExpand = TRUE
    )
  }
}

saveWorkbook(wb, out_file, overwrite = TRUE)

cat("\n================ SUPPLEMENTARY TABLE S1 COMPLETE ================\n")
cat("Output: ", out_file, "\n", sep = "")
cat("Subtables: S1A-S1H\n")
cat("S1A global tests: ", nrow(s1a), " rows\n", sep = "")
cat("S1B BeWo contrasts: ", nrow(s1b), " rows\n", sep = "")
cat("S1C exposure contrasts: ", nrow(s1c), " rows\n", sep = "")
cat("S1D interaction tests: ", nrow(s1d), " rows\n", sep = "")
cat("S1E interaction contrasts: ", nrow(s1e), " rows\n", sep = "")
cat("S1F all-pairwise cell-line contrasts: ", nrow(s1f), " rows\n", sep = "")
cat("S1G cell-line summaries: ", nrow(s1g), " rows\n", sep = "")
cat("S1H condition summaries: ", nrow(s1h), " rows\n", sep = "")
cat("=================================================================\n")
