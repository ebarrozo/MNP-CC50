# ==============================================================================
# make_Table1_v6_3_FINAL.R
# Reproducibly generate consolidated manuscript Table 1 from final v6.3 outputs
#
# Input:
#   v6_3_final_BeWo_outputs/V6_3_CC50_ByCondition_Descriptive.csv
#
# Outputs:
#   manuscript_outputs/Table1_Final_Consolidated_CC50.csv
#   manuscript_outputs/Table1_Final_Consolidated_CC50.docx
#
# Table:
#   Median boundary-assigned CC50 (ug/mL)
#   Explicit censoring notation: <=0.2 and >20
#   Polymer grouped by rows; timepoint x size across columns
# ==============================================================================

# Run this script from the repository root.
project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(project_dir, "scripts")) ||
    !dir.exists(file.path(project_dir, "metadata"))) {
  stop("Run this script from the repository root (the directory containing scripts/ and metadata/).")
}


suppressPackageStartupMessages({
  library(tidyverse)
  library(flextable)
  library(officer)
})

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

input_file <- file.path(
  project_dir,
  "results", "v6_3",
  "V6_3_CC50_ByCondition_Descriptive.csv"
)

out_dir <- file.path(project_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

csv_out <- file.path(out_dir, "Table1_Final_Consolidated_CC50.csv")
docx_out <- file.path(out_dir, "Table1_Final_Consolidated_CC50.docx")

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

# ------------------------------------------------------------------------------
# 1. Load final condition-level descriptive results
# ------------------------------------------------------------------------------

dat <- readr::read_csv(input_file, show_col_types = FALSE)

required <- c(
  "cell_line", "Treatment", "Size", "Timepoint",
  "n", "median_boundary", "n_lower", "n_upper"
)

missing_cols <- setdiff(required, names(dat))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 2. Ordering and display labels
# ------------------------------------------------------------------------------

polymer_levels <- c("Mix", "PE", "PS", "PVC")
cell_levels <- c("BeWo", "HTR8", "JEG3", "THP-1", "Jurkat")
time_levels <- c(3, 6, 24, 48)
size_levels <- c("Small", "Large")

dat <- dat %>%
  mutate(
    Treatment = factor(Treatment, levels = polymer_levels),
    cell_line = factor(cell_line, levels = cell_levels),
    Timepoint = as.numeric(as.character(Timepoint)),
    Size = factor(Size, levels = size_levels),

    Cell_line_display = recode(
      as.character(cell_line),
      "HTR8" = "HTR-8/SVneo",
      "JEG3" = "JEG-3",
      .default = as.character(cell_line)
    )
  )

# ------------------------------------------------------------------------------
# 3. Format median CC50 with explicit assay-boundary censoring
# ------------------------------------------------------------------------------

format_cc50 <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    abs(x - 0.2) < 1e-10 ~ "\u22640.2",
    abs(x - 20) < 1e-10 ~ ">20",
    TRUE ~ sprintf("%.1f", x)
  )
}

dat <- dat %>%
  mutate(
    CC50_display = format_cc50(median_boundary),
    Column = paste0(Timepoint, " h ", as.character(Size))
  )

# ------------------------------------------------------------------------------
# 4. Build one consolidated table
# ------------------------------------------------------------------------------

table1 <- dat %>%
  select(
    Polymer = Treatment,
    `Cell line` = Cell_line_display,
    Timepoint,
    Size,
    CC50_display
  ) %>%
  arrange(
    factor(Polymer, levels = polymer_levels),
    factor(`Cell line`,
           levels = c("BeWo", "HTR-8/SVneo", "JEG-3", "THP-1", "Jurkat")),
    Timepoint,
    Size
  ) %>%
  mutate(
    Column = paste0(Timepoint, " h ", Size)
  ) %>%
  select(Polymer, `Cell line`, Column, CC50_display) %>%
  pivot_wider(
    names_from = Column,
    values_from = CC50_display
  ) %>%
  select(
    Polymer,
    `Cell line`,
    `3 h Small`, `3 h Large`,
    `6 h Small`, `6 h Large`,
    `24 h Small`, `24 h Large`,
    `48 h Small`, `48 h Large`
  )

# Save machine-readable manuscript table.
readr::write_csv(table1, csv_out)

# ------------------------------------------------------------------------------
# 5. Replicate-count audit for footnote
# ------------------------------------------------------------------------------

n_audit <- dat %>%
  distinct(cell_line, Timepoint, n) %>%
  arrange(cell_line, Timepoint)

readr::write_csv(
  n_audit,
  file.path(out_dir, "Table1_BiologicalReplicate_Audit.csv")
)

# ------------------------------------------------------------------------------
# 6. Build manuscript-ready flextable
# ------------------------------------------------------------------------------

ft <- flextable(table1)

# Two-level column headers: Timepoint with Small/Large beneath.
header_map <- data.frame(
  key = names(table1),
  line1 = c(
    "Polymer", "Cell line",
    "3 h", "3 h",
    "6 h", "6 h",
    "24 h", "24 h",
    "48 h", "48 h"
  ),
  line2 = c(
    "", "",
    "Small", "Large",
    "Small", "Large",
    "Small", "Large",
    "Small", "Large"
  ),
  stringsAsFactors = FALSE
)

ft <- set_header_df(
  ft,
  mapping = header_map,
  key = "key"
)

ft <- merge_h(ft, part = "header")
ft <- merge_v(ft, j = 1, part = "body")

# General formatting
ft <- theme_booktabs(ft)
ft <- font(ft, fontname = "Arial", part = "all")
ft <- fontsize(ft, size = 8.5, part = "body")
ft <- fontsize(ft, size = 8.5, part = "header")
ft <- bold(ft, part = "header")
ft <- bold(ft, j = 1, part = "body")

# Alignment
ft <- align(ft, j = 1, align = "center", part = "all")
ft <- align(ft, j = 2, align = "left", part = "all")
ft <- align(ft, j = 3:10, align = "center", part = "all")
ft <- valign(ft, valign = "center", part = "all")

# Compact dimensions
ft <- width(ft, j = 1, width = 0.65)
ft <- width(ft, j = 2, width = 1.15)
ft <- width(ft, j = 3:10, width = 0.78)

ft <- padding(
  ft,
  padding.top = 2,
  padding.bottom = 2,
  padding.left = 3,
  padding.right = 3,
  part = "all"
)

# Repeat headers on page break if journal/editor opens it outside landscape.
ft <- set_table_properties(
  ft,
  layout = "fixed",
  width = 1
)

# ------------------------------------------------------------------------------
# 7. Title and footnote
# ------------------------------------------------------------------------------

title_text <- paste(
  "Table 1. Median CC50 values across placental and immune cell models",
  "following micro- and nanoplastic exposure"
)

footnote_text <- paste0(
  "Values are median boundary-assigned CC50 estimates (\u00b5g/mL) across ",
  "QC-valid independent biological experiments. Values at the limits of the ",
  "tested concentration range are reported as \u22640.2 or >20 \u00b5g/mL. ",
  "BeWo and Jurkat include n=3 biological replicates for all conditions; ",
  "HTR-8/SVneo and THP-1 include n=2 after plate-level QC; JEG-3 includes ",
  "n=2 at 3 and 6 h and n=3 at 24 and 48 h. Small and Large denote the ",
  "polymer-specific particle-size classes defined in Methods. ",
  "CC50, 50% cytotoxic concentration."
)

ft <- set_caption(
  ft,
  caption = as_paragraph(as_b(title_text)),
  autonum = NULL
)

ft <- add_footer_lines(
  ft,
  values = footnote_text
)

ft <- fontsize(ft, size = 8, part = "footer")
ft <- align(ft, align = "left", part = "footer")

# ------------------------------------------------------------------------------
# 8. Export Word document
# ------------------------------------------------------------------------------

doc <- read_docx()

doc <- body_end_section_landscape(doc)

doc <- body_add_flextable(
  x = doc,
  value = ft,
  align = "center"
)

print(doc, target = docx_out)

# ------------------------------------------------------------------------------
# 9. Console audit
# ------------------------------------------------------------------------------

cat("\n==================== TABLE 1 COMPLETE ====================\n")
cat("Input: ", input_file, "\n", sep = "")
cat("CSV: ", csv_out, "\n", sep = "")
cat("Word: ", docx_out, "\n", sep = "")
cat("Rows: ", nrow(table1), " (expected 20)\n", sep = "")
cat("Condition cells: ", nrow(table1) * 8, " (expected 160)\n", sep = "")
cat("Summary: median boundary-assigned CC50; no SEM/range\n")
cat("Censor notation: <=0.2 and >20 ug/mL\n")
cat("==========================================================\n")
