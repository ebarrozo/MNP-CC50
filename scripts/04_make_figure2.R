# ==============================================================================
# make_Figure2_v6_3_FINAL_v3.R
# Final Figure 2
#
# Panel A:
#   Condition-level median boundary-assigned CC50 heatmap
#   Facets = particle size (rows) x polymer (columns)
#   X-axis = exposure duration
#
# Panel B:
#   Proportion of replicate-level CC50 estimates classified as
#   <=0.2, interpolated 0.2-20, or >20 ug/mL
#
# Inputs:
#   v6_3_final_BeWo_outputs/V6_3_CC50_ByCondition_Descriptive.csv
#   v5_2_outputs/CC50_ReplicateValues.csv
#
# Outputs:
#   manuscript_outputs/Figure2_CC50_Landscape_FINAL.pdf
#   manuscript_outputs/Figure2_CC50_Landscape_FINAL.png
#   manuscript_outputs/Figure2_CC50_Landscape_FINAL.tiff
#   manuscript_outputs/Figure2A_PlotData.csv
#   manuscript_outputs/Figure2B_PlotData.csv
# ==============================================================================

# Run this script from the repository root.
project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(project_dir, "scripts")) ||
    !dir.exists(file.path(project_dir, "metadata"))) {
  stop("Run this script from the repository root (the directory containing scripts/ and metadata/).")
}


suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------------------------
# 0. Paths and constants
# ------------------------------------------------------------------------------

condition_file <- file.path(
  project_dir,
  "results", "v6_3",
  "V6_3_CC50_ByCondition_Descriptive.csv"
)

replicate_file <- file.path(
  project_dir,
  "results", "v5_2",
  "CC50_ReplicateValues.csv"
)

out_dir <- file.path(project_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(condition_file, replicate_file)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

LOW_CAP <- 0.2
HIGH_CAP <- 20

# Individual polymers first, then mixture.
polymer_levels <- c("PE", "PS", "PVC", "Mix")

# Preserve biologically intuitive grouping:
# trophoblast models followed by immune-cell models.
cell_display_levels <- c(
  "BeWo",
  "HTR-8/SVneo",
  "JEG-3",
  "THP-1",
  "Jurkat"
)

time_levels <- c(3, 6, 24, 48)
size_levels <- c("Small", "Large")

# ------------------------------------------------------------------------------
# 1. Panel A: median CC50 landscape
# ------------------------------------------------------------------------------

cond <- readr::read_csv(condition_file, show_col_types = FALSE)

required_cond <- c(
  "cell_line", "Treatment", "Size", "Timepoint",
  "n", "median_boundary", "n_lower", "n_interpolated", "n_upper"
)

missing_cond <- setdiff(required_cond, names(cond))
if (length(missing_cond) > 0) {
  stop(
    "Missing condition-summary columns: ",
    paste(missing_cond, collapse = ", ")
  )
}

cond <- cond %>%
  mutate(
    cell_line = recode(
      as.character(cell_line),
      "HTR8" = "HTR-8/SVneo",
      "JEG3" = "JEG-3",
      .default = as.character(cell_line)
    ),
    # Reverse here because ggplot draws first factor level at bottom.
    cell_line = factor(cell_line, levels = rev(cell_display_levels)),
    Treatment = factor(Treatment, levels = polymer_levels),
    Timepoint = as.numeric(as.character(Timepoint)),
    Timepoint_f = factor(Timepoint, levels = time_levels),
    Size = factor(Size, levels = size_levels)
  )

pA <- ggplot(
  cond,
  aes(
    x = Timepoint_f,
    y = cell_line,
    fill = median_boundary
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  facet_grid(
    rows = vars(Size),
    cols = vars(Treatment),
    switch = "y"
  ) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(LOW_CAP, HIGH_CAP),
    breaks = c(0.2, 5, 10, 15, 20),
    labels = c("≤0.2", "5", "10", "15", ">20"),
    name = expression("Median CC"[50]*"\n(µg/mL)")
  ) +
  scale_x_discrete(
    labels = paste0(time_levels, " h")
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "A",
    subtitle = expression("CC"[50]*" landscape by exposure condition")
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0
    ),
    plot.subtitle = element_text(
      face = "bold",
      size = 10.5
    ),
    axis.text.x = element_text(size = 8.5),
    axis.text.y = element_text(size = 9),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      face = "bold",
      size = 10
    ),
    strip.text.y.left = element_text(
      face = "bold",
      angle = 0,
      size = 9
    ),
    strip.placement = "outside",
    panel.spacing = unit(0.7, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )

# ------------------------------------------------------------------------------
# 2. Panel B: replicate-level assay-range distribution
# ------------------------------------------------------------------------------

reps <- readr::read_csv(replicate_file, show_col_types = FALSE)

if (!"cc50_censor" %in% names(reps) &&
    "CC50_Censor" %in% names(reps)) {
  reps <- reps %>%
    rename(cc50_censor = CC50_Censor)
}

required_rep <- c(
  "cell_line",
  "cc50",
  "cc50_censor"
)

missing_rep <- setdiff(required_rep, names(reps))
if (length(missing_rep) > 0) {
  stop(
    "Missing replicate-level columns: ",
    paste(missing_rep, collapse = ", ")
  )
}

reps <- reps %>%
  mutate(
    cell_line = recode(
      as.character(cell_line),
      "HTR8" = "HTR-8/SVneo",
      "JEG3" = "JEG-3",
      .default = as.character(cell_line)
    ),
    cell_line = factor(
      cell_line,
      levels = cell_display_levels
    ),
    cc50 = suppressWarnings(as.numeric(cc50)),

    CC50_Category = case_when(
      str_detect(
        toupper(as.character(cc50_censor)),
        "LOWER|<=|BELOW"
      ) ~ "≤0.2",

      str_detect(
        toupper(as.character(cc50_censor)),
        "UPPER|>|ABOVE"
      ) ~ ">20",

      str_detect(
        toupper(as.character(cc50_censor)),
        "INTERP"
      ) ~ "0.2–20",

      is.finite(cc50) &
        abs(cc50 - LOW_CAP) < 1e-10 ~ "≤0.2",

      is.finite(cc50) &
        abs(cc50 - HIGH_CAP) < 1e-10 ~ ">20",

      is.finite(cc50) ~ "0.2–20",

      TRUE ~ NA_character_
    ),

    # Factor order chosen so stacked bar reads:
    # low CC50 at bottom -> within range -> upper-censored at top.
    CC50_Category = factor(
      CC50_Category,
      levels = c(">20", "0.2–20", "≤0.2")
    )
  ) %>%
  filter(
    is.finite(cc50),
    !is.na(CC50_Category)
  )

censor <- reps %>%
  count(
    cell_line,
    CC50_Category,
    name = "n"
  ) %>%
  group_by(cell_line) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup()

pB <- ggplot(
  censor,
  aes(
    x = cell_line,
    y = proportion,
    fill = CC50_Category
  )
) +
  geom_col(width = 0.70) +
  geom_text(
    aes(
      label = if_else(
        proportion >= 0.10,
        percent(proportion, accuracy = 1),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    size = 3.0,
    fontface = "bold",
    color = "white"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_viridis_d(
    option = "D",
    begin = 0.15,
    end = 0.85,
    direction = -1,
    breaks = c("≤0.2", "0.2–20", ">20"),
    name = expression("Replicate CC"[50])
  ) +
  labs(
    x = NULL,
    y = expression(
      "Proportion of replicate-level CC"[50]*" estimates"
    ),
    title = "B",
    subtitle = "Distribution relative to assay range"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0
    ),
    plot.subtitle = element_text(
      face = "bold",
      size = 10.5
    ),
    axis.text.x = element_text(
      size = 9,
      angle = 25,
      hjust = 1
    ),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )

# ------------------------------------------------------------------------------
# 3. Export exact plotting data
# ------------------------------------------------------------------------------

readr::write_csv(
  cond %>%
    select(
      cell_line,
      Treatment,
      Size,
      Timepoint,
      n,
      median_boundary,
      n_lower,
      n_interpolated,
      n_upper
    ),
  file.path(
    out_dir,
    "Figure2A_PlotData.csv"
  )
)

readr::write_csv(
  censor,
  file.path(
    out_dir,
    "Figure2B_PlotData.csv"
  )
)

# ------------------------------------------------------------------------------
# 4. Assemble final Figure 2
# ------------------------------------------------------------------------------

fig2 <- pA / pB +
  plot_layout(
    heights = c(1.6, 1)
  )

ggsave(
  file.path(
    out_dir,
    "Figure2_CC50_Landscape_FINAL.pdf"
  ),
  fig2,
  width = 10.0,
  height = 8.7,
  units = "in",
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Figure2_CC50_Landscape_FINAL.png"
  ),
  fig2,
  width = 10.0,
  height = 8.7,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Figure2_CC50_Landscape_FINAL.tiff"
  ),
  fig2,
  width = 10.0,
  height = 8.7,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# Separate panel exports.
ggsave(
  file.path(
    out_dir,
    "Figure2A_CC50_Heatmap_FINAL.pdf"
  ),
  pA,
  width = 10.0,
  height = 5.3,
  units = "in",
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Figure2B_CC50_AssayRange_FINAL.pdf"
  ),
  pB,
  width = 7.0,
  height = 4.0,
  units = "in",
  bg = "white"
)

# ------------------------------------------------------------------------------
# 5. Console audit
# ------------------------------------------------------------------------------

cat("\n==================== FIGURE 2 FINAL ====================\n")
cat(
  "Heatmap cells: ",
  nrow(cond),
  " / 160 expected\n",
  sep = ""
)
cat(
  "Replicate CC50 observations: ",
  nrow(reps),
  " / 400 expected\n",
  sep = ""
)
cat("Polymer order: PE, PS, PVC, Mix\n")
cat("Boundary triangles: removed\n")
cat(
  "Panel B y-axis: proportion of replicate-level CC50 estimates\n"
)
cat(
  "Output: manuscript_outputs/Figure2_CC50_Landscape_FINAL.*\n"
)
cat("========================================================\n")
