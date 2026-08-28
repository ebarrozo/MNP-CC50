# MNP CC50 cytotoxicity analysis

Reproducibility repository for the manuscript **“Polymer identity and cell type shape microplastic and nanoplastic cytotoxicity across human placental and immune cell lines.”**

This workflow analyzes Ethidium homodimer-1 (EthD-1) fluorescence after exposure of five human placental and immune cell lines to fluorescent polystyrene (PS), polyethylene (PE), polyvinyl chloride (PVC), or a 1:1:1 PS:PE:PVC mixture.

## Experimental design

- Cell lines: BeWo, HTR-8/SVneo, JEG-3, THP-1, Jurkat
- Polymer treatments: PS, PE, PVC, Mix (1:1:1 PS:PE:PVC by mass)
- Nominal particle-size classes:
  - Small: 50–100 nm, depending on polymer
  - Large: 10–20 µm, depending on polymer
- Final well concentrations: **0.2, 2, and 20 µg/mL**
- Exposure durations: 3, 6, 24, and 48 h
- Three technical wells per dose within each independent experiment
- Nominally three independent experiments per cell line/timepoint; plate-level QC determines the analysis-ready replicate count

### Important dose note

The historical plate-layout workbook retains the original working-suspension labels `1`, `10`, and `100`. Because 50 µL of each working suspension was added to 200 µL already present in the well, the analytically correct **final well concentrations are 0.2, 2, and 20 µg/mL**. All analysis scripts use the corrected final concentrations.

## CC50 definition

CC50 is the assay-defined concentration associated with 50% EthD-1-defined membrane-integrity loss relative to plate-specific media and methanol controls.

Technical triplicates are averaged before CC50 estimation. If adjacent tested concentrations bracket 50% cell death, CC50 is estimated by linear interpolation. If 50% is not bracketed:

- all tested concentrations at or above 50% death → true CC50 `≤0.2 µg/mL`, represented numerically as `0.2` in the primary boundary-assigned model;
- all tested concentrations below 50% death → true CC50 `>20 µg/mL`, represented numerically as `20` in the primary boundary-assigned model.

Because many estimates are censored at the assay boundaries, the primary linear mixed-effects model is accompanied by a censoring-aware ordinal mixed-effects sensitivity analysis.

## Repository structure

```text
.
├── README.md
├── scripts/
│   ├── 01_preprocess_plate_data.R
│   ├── 02_fit_cc50_models.R
│   ├── 03_make_table1.R
│   ├── 04_make_figure2.R
│   ├── 05_make_figure3.R
│   ├── 06_make_supplementary_table_s1.R
│   ├── 07_particle_interference_diagnostics.R
│   ├── 08_particle_interference_condition_models.R
│   └── 09_make_supplementary_table_s2.R
├── metadata/
│   ├── key_plate_layout.xlsx
│   └── README.md
├── data/
│   ├── README.md
│   ├── raw/
│   └── no_cell/
│       └── 08.26.26-EthD1_mnponly.xlsx
├── results/                 # generated intermediate/model outputs
├── outputs/                 # generated manuscript tables and figures
└── docs/
    └── analysis_workflow.md
```

## Required R packages

Core workflow:

```r
install.packages(c(
  "tidyverse", "readxl", "stringr",
  "lme4", "lmerTest", "emmeans", "ordinal",
  "broom", "patchwork", "scales",
  "openxlsx", "flextable", "officer"
))
```

For the HC3 robust-standard-error sensitivity analysis in script 08:

```r
install.packages(c("sandwich", "lmtest"))
```

The scripts should be run from the **repository root**. No user-specific absolute paths are required.

## Running the analysis

Place the 60 raw BMG MARS `.xlsx` exports in `data/raw/` using their original filenames, then run from a fresh R session:

```r
source("scripts/01_preprocess_plate_data.R")
source("scripts/02_fit_cc50_models.R")
source("scripts/03_make_table1.R")
source("scripts/04_make_figure2.R")
source("scripts/05_make_figure3.R")
source("scripts/06_make_supplementary_table_s1.R")
source("scripts/07_particle_interference_diagnostics.R")
source("scripts/08_particle_interference_condition_models.R")
source("scripts/09_make_supplementary_table_s2.R")
```

See [`docs/analysis_workflow.md`](docs/analysis_workflow.md) for analysis details and expected outputs.

## Plate-level quality control

The primary normalization is

```text
100 × (sample − media) / (methanol − media)
```

Primary plate validity requires methanol fluorescence to exceed media fluorescence. Control separation is additionally characterized using the strictly standardized mean difference (SSMD). Plates with SSMD <1 are retained but flagged in the primary workflow; a stricter sensitivity analysis includes only plates with SSMD ≥1. Technical wells are not treated as independent statistical observations.

## Statistical analysis

The primary model is an additive linear mixed-effects model of boundary-assigned replicate-level CC50:

```text
CC50 ~ cell line + polymer treatment + particle-size class + exposure duration + (1 | biological replicate)
```

BeWo is the reference cell line. Secondary models evaluate one cell-line interaction at a time. P values are adjusted using the Benjamini-Hochberg false-discovery-rate procedure within prespecified inferential families.

The censoring-aware sensitivity analysis classifies CC50 as `≤0.2`, `0.2–20`, or `>20 µg/mL` and fits a cumulative-link mixed model using the same additive fixed effects and biological-replicate random intercept.

## Particle-associated fluorescence control

The fluorescent MNP preparations generated measurable signal within the EthD-1 acquisition window. The no-cell experiment is therefore used to characterize potential optical interference.

The no-cell controls were measured **without the wash steps used before EthD-1 acquisition in cell-containing wells**. Accordingly, they are treated as a conservative characterization of potential particle-associated fluorescence rather than a direct estimate of residual post-wash particle signal. No-cell fluorescence is **not subtracted** from the primary assay because its contribution in cell-containing wells cannot be assumed to be linearly additive.

Scripts 07 and 08 instead test whether particle-associated fluorescence systematically tracks the original apparent cell-death signal across the 24 unique polymer × particle-size × concentration conditions. The public workflow intentionally excludes the exploratory direct-subtraction analysis because that additivity assumption was not supported.

## Data availability

The no-cell optical-interference control and plate-layout metadata are included here. The complete set of raw BMG MARS plate-reader exports must also be placed in `data/raw/` for end-to-end reproduction of the primary analysis.

If the raw files are deposited in an external repository rather than GitHub, add the accession/DOI here before publication:

> **Raw-data accession:** [TO ADD]

## Code provenance

The numbered public scripts correspond to the finalized manuscript analysis. Version-like strings retained in some generated filenames (for example, `V6_3_...`) are provenance labels from the locked analysis and do not indicate that earlier exploratory workflows must be run.

## License

A software/data license has not yet been selected. Add the appropriate license before making the repository public.
