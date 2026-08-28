# ==============================================================================
# make_Figure3_v6_3_FINAL_v4.R
# Final inferential summary figure for the corrected MNP CC50 analysis
#
# Panel A: Global adjusted effects in the primary additive LMM and
#          censoring-aware ordinal CLMM.
# Panel B: Primary LMM cell-line contrasts versus BeWo.
# Panel C: Ordinal CLMM cell-line contrasts versus BeWo.
# Panel D: Mixture-versus-individual-polymer contrasts in both models.
#
# All inputs are FINAL v6.3 exported results. No models are refit here.
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
# 0. Paths
# ------------------------------------------------------------------------------

res_dir <- file.path(project_dir, "results", "v6_3")
out_dir <- file.path(project_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files_needed <- c(
  "V6_3_PRIMARY_Additive_LMM_TypeIII.csv",
  "V6_3_SENS_OrdinalCLMM_LRT_Tests.csv",
  "V6_3_PRIMARY_CellLine_vs_BeWo.csv",
  "V6_3_SENS_OrdinalCLMM_CellLine_vs_BeWo.csv",
  "V6_3_PRIMARY_Treatment_Pairwise.csv",
  "V6_3_SENS_OrdinalCLMM_Treatment_Pairwise.csv"
)

for (f in files_needed) {
  fp <- file.path(res_dir, f)
  if (!file.exists(fp)) stop("Missing required v6.3 output: ", fp)
}

# ------------------------------------------------------------------------------
# 1. Panel A: omnibus/global adjusted effects
# ------------------------------------------------------------------------------

global_lmm <- readr::read_csv(
  file.path(res_dir, "V6_3_PRIMARY_Additive_LMM_TypeIII.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    Term,
    Model = "Primary LMM",
    q_BH
  )

global_ord <- readr::read_csv(
  file.path(res_dir, "V6_3_SENS_OrdinalCLMM_LRT_Tests.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    Term,
    Model = "Ordinal sensitivity",
    q_BH
  )

global <- bind_rows(global_lmm, global_ord) %>%
  mutate(
    Term = recode(
      Term,
      "cell_line" = "Cell line",
      "Treatment" = "Polymer",
      "Size" = "Particle size",
      "Timepoint" = "Exposure duration"
    ),
    Term = factor(
      Term,
      levels = rev(c(
        "Cell line",
        "Polymer",
        "Particle size",
        "Exposure duration"
      ))
    ),
    Model = factor(
      Model,
      levels = c("Primary LMM", "Ordinal sensitivity")
    ),
    minuslog10q = -log10(q_BH)
  )

pA <- ggplot(
  global,
  aes(x = minuslog10q, y = Term, shape = Model)
) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = 2,
    linewidth = 0.55
  ) +
  geom_point(size = 3.2) +
  scale_x_continuous(
    limits = c(0, max(global$minuslog10q) + 1.0),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = expression(-log[10]("BH-adjusted q")),
    y = NULL,
    shape = NULL,
    title = "A",
    subtitle = "Global adjusted effects"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(face = "bold", size = 10.5),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    legend.text = element_text(size = 8.5),
    plot.margin = margin(5, 8, 5, 5)
  )

# ------------------------------------------------------------------------------
# 2. Panel B: primary LMM cell-line contrasts versus BeWo
# ------------------------------------------------------------------------------

cell_lmm <- readr::read_csv(
  file.path(res_dir, "V6_3_PRIMARY_CellLine_vs_BeWo.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    test_cell = str_replace_all(test_cell, fixed("(THP-1)"), "THP-1"),
    test_cell = recode(
      test_cell,
      "HTR8" = "HTR-8/SVneo",
      "JEG3" = "JEG-3"
    ),
    test_cell = factor(
      test_cell,
      levels = rev(c("HTR-8/SVneo", "JEG-3", "THP-1", "Jurkat"))
    ),
    sig = q_BH < 0.05,
    q_label = case_when(
      q_BH < 0.001 ~ "q<0.001",
      TRUE ~ paste0("q=", formatC(q_BH, format = "f", digits = 3))
    )
  )

pB <- ggplot(
  cell_lmm,
  aes(x = Estimate, y = test_cell)
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.55
  ) +
  geom_errorbar(
    aes(xmin = lowerCL, xmax = upperCL),
    width = 0.15,
    linewidth = 0.75
  ) +
  geom_point(size = 3.2) +
  labs(
    x = expression(Delta*"CC"[50]*" vs BeWo (µg/mL)"),
    y = NULL,
    title = "B",
    subtitle = "Primary additive LMM"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(face = "bold", size = 10.5),
    axis.text = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5)
  )

# ------------------------------------------------------------------------------
# 3. Panel C: ordinal CLMM cell-line contrasts versus BeWo
# ------------------------------------------------------------------------------

cell_ord <- readr::read_csv(
  file.path(res_dir, "V6_3_SENS_OrdinalCLMM_CellLine_vs_BeWo.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    test_cell = str_replace_all(test_cell, fixed("(THP-1)"), "THP-1"),
    test_cell = recode(
      test_cell,
      "HTR8" = "HTR-8/SVneo",
      "JEG3" = "JEG-3"
    ),
    test_cell = factor(
      test_cell,
      levels = rev(c("HTR-8/SVneo", "JEG-3", "THP-1", "Jurkat"))
    ),
    sig = q_BH < 0.05
  )

pC <- ggplot(
  cell_ord,
  aes(x = OR_higher_CC50_category, y = test_cell)
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.55
  ) +
  geom_errorbar(
    aes(xmin = OR_lower95, xmax = OR_upper95),
    width = 0.15,
    linewidth = 0.75
  ) +
  geom_point(size = 3.2) +
  scale_x_log10(
    breaks = c(0.01, 0.1, 1, 10, 100, 300),
    labels = c("0.01", "0.1", "1", "10", "100", "300")
  ) +
  labs(
    x = expression("OR for higher CC"[50]*" category"),
    y = NULL,
    title = "C",
    subtitle = "Censoring-aware ordinal sensitivity"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(face = "bold", size = 10.5),
    axis.text = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5)
  )

# ------------------------------------------------------------------------------
# 4. Panel D: mixture versus individual polymers
# ------------------------------------------------------------------------------

trt_lmm <- readr::read_csv(
  file.path(res_dir, "V6_3_PRIMARY_Treatment_Pairwise.csv"),
  show_col_types = FALSE
) %>%
  filter(contrast %in% c("Mix - PE", "Mix - PS", "Mix - PVC")) %>%
  transmute(
    contrast,
    Comparator = str_remove(contrast, "^Mix - "),
    Model = "Primary LMM",
    estimate = Estimate,
    lower = lowerCL,
    upper = upperCL,
    q_BH,
    metric = "Difference"
  )

trt_ord <- readr::read_csv(
  file.path(res_dir, "V6_3_SENS_OrdinalCLMM_Treatment_Pairwise.csv"),
  show_col_types = FALSE
) %>%
  filter(contrast %in% c("Mix - PE", "Mix - PS", "Mix - PVC")) %>%
  transmute(
    contrast,
    Comparator = str_remove(contrast, "^Mix - "),
    Model = "Ordinal sensitivity",
    estimate = OR_higher_CC50_category,
    lower = OR_lower95,
    upper = OR_upper95,
    q_BH,
    metric = "Odds ratio"
  )

trt_lmm$Comparator <- factor(trt_lmm$Comparator, levels = rev(c("PE", "PS", "PVC")))
trt_ord$Comparator <- factor(trt_ord$Comparator, levels = rev(c("PE", "PS", "PVC")))

pD1 <- ggplot(
  trt_lmm,
  aes(x = estimate, y = Comparator)
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.55
  ) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0.15,
    linewidth = 0.75
  ) +
  geom_point(size = 3.1) +
  labs(
    x = expression(Delta*"CC"[50]*": Mix minus polymer (µg/mL)"),
    y = NULL,
    title = "D",
    subtitle = "Primary LMM"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(face = "bold", size = 10.5),
    axis.text = element_text(size = 8.5),
    plot.margin = margin(0, 5, 0, 5)
  )

pD2 <- ggplot(
  trt_ord,
  aes(x = estimate, y = Comparator)
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.55
  ) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0.15,
    linewidth = 0.75
  ) +
  geom_point(size = 3.1) +
  scale_x_log10(
    breaks = c(1, 2, 4, 8),
    labels = c("1", "2", "4", "8")
  ) +
  labs(
    x = expression("OR for higher CC"[50]*" category"),
    y = NULL,
    subtitle = "Ordinal sensitivity"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    axis.text = element_text(size = 8.5),
    plot.margin = margin(0, 5, 0, 5)
  )

pD <- pD1 | pD2

# ------------------------------------------------------------------------------
# 5. Export plotted data for reproducibility
# ------------------------------------------------------------------------------

readr::write_csv(global, file.path(out_dir, "Figure3A_GlobalEffects_PlotData.csv"))
readr::write_csv(cell_lmm, file.path(out_dir, "Figure3B_Primary_CellLine_PlotData.csv"))
readr::write_csv(cell_ord, file.path(out_dir, "Figure3C_Ordinal_CellLine_PlotData.csv"))
readr::write_csv(trt_lmm, file.path(out_dir, "Figure3D_Primary_Treatment_PlotData.csv"))
readr::write_csv(trt_ord, file.path(out_dir, "Figure3D_Ordinal_Treatment_PlotData.csv"))

# ------------------------------------------------------------------------------
# 6. Assemble/export final figure
# ------------------------------------------------------------------------------

top <- pA | pB
bottom <- pC | pD

fig3 <- top / bottom +
  plot_layout(heights = c(1, 1.05))

ggsave(
  file.path(out_dir, "Figure3_InferentialSummary_FINAL.pdf"),
  fig3,
  width = 11.0,
  height = 8.0,
  units = "in",
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure3_InferentialSummary_FINAL.png"),
  fig3,
  width = 11.0,
  height = 8.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure3_InferentialSummary_FINAL.tiff"),
  fig3,
  width = 11.0,
  height = 8.0,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# Separate vector panels
ggsave(file.path(out_dir, "Figure3A_GlobalEffects_FINAL.pdf"),
       pA, width = 5.2, height = 3.8, units = "in", bg = "white")
ggsave(file.path(out_dir, "Figure3B_Primary_CellLine_FINAL.pdf"),
       pB, width = 5.2, height = 3.8, units = "in", bg = "white")
ggsave(file.path(out_dir, "Figure3C_Ordinal_CellLine_FINAL.pdf"),
       pC, width = 5.2, height = 3.8, units = "in", bg = "white")
ggsave(file.path(out_dir, "Figure3D_Treatment_FINAL.pdf"),
       pD, width = 5.2, height = 3.8, units = "in", bg = "white")

# ------------------------------------------------------------------------------
# 7. Console audit
# ------------------------------------------------------------------------------

cat("\n==================== FIGURE 3 FINAL ====================\n")
cat("Panel A global effects: ", nrow(global), " points (8 expected)\n", sep = "")
cat("Panel B BeWo primary contrasts: ", nrow(cell_lmm), " (4 expected)\n", sep = "")
cat("Panel C BeWo ordinal contrasts: ", nrow(cell_ord), " (4 expected)\n", sep = "")
cat("Panel D Mix primary contrasts: ", nrow(trt_lmm), " (3 expected)\n", sep = "")
cat("Panel D Mix ordinal contrasts: ", nrow(trt_ord), " (3 expected)\n", sep = "")
cat("No secondary interaction contrasts are plotted.\n")
cat("Panel A q-value text labels and horizontal stems: removed for readability.\n")
cat("Panels B/C/D points: solid; significance conveyed by CI/q-values in text.\n")
cat("Panel D subplots labeled Primary LMM and Ordinal sensitivity; shared descriptive title moved to legend.\n")
cat("Output: manuscript_outputs/Figure3_InferentialSummary_FINAL.*\n")
cat("========================================================\n")
