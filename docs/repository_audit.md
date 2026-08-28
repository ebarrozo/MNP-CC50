# Repository audit and cleanup

This repository was prepared from the manuscript analysis scripts supplied for public release.

## Critical issues found in the source folder

### 1. User-specific absolute paths

All nine source scripts contained a hard-coded `setwd()` path to a local Box directory. These were removed. Public scripts now use paths relative to the repository root and explicitly stop if run from the wrong directory.

### 2. Hidden dependency on an excluded exploratory v6.4 analysis

The original particle-interference scripts depended on:

```text
v6_4_particle_fluorescence_sensitivity_outputs/
V6_4_NoCell_ParticleFluorescence_ByCondition.csv
```

The v6.4 workflow was the exploratory direct-subtraction analysis that was ultimately not used in the manuscript and was not included in the supplied GitHub folder. This made the original v6.5/v6.6/S2 workflow non-reproducible.

The public script `07_particle_interference_diagnostics.R` now reads the original no-cell workbook directly and generates the required 24-condition fluorescence summary. The direct-subtraction workflow is not required or exposed.

### 3. Repeated-observation interference LMM

The original v6.5 diagnostic included an LMM that repeated the same 24 particle-only fluorescence values across many cell-containing observations, producing a much larger nominal denominator degrees of freedom than the no-cell experiment supports. This model is not used in the manuscript's final inference and has been removed from the public workflow. The public analysis emphasizes the 24-condition correlations and the final condition-level regressions.

### 4. Script naming/version clutter

Public script names were standardized to the numbered workflow `01` through `09`. The S1 script no longer exposes `PATCHED` in its filename. Internal `V6_3_...` output names are retained as locked-analysis provenance labels.

### 5. Historical dose labels

The plate key contains historical 1/10/100 working-suspension labels. The public README and metadata documentation explicitly explain that the final well concentrations were 0.2/2/20 µg/mL and that the analysis uses those corrected values.

### 6. macOS archive artifacts

`__MACOSX` and `._*` files from the uploaded ZIP were excluded from the cleaned repository.

## Files added

- Root `README.md`
- `docs/analysis_workflow.md`
- `docs/repository_audit.md`
- `metadata/README.md`
- `data/README.md`
- `results/README.md`
- `outputs/README.md`
- `.gitignore`
- no-cell source workbook under `data/no_cell/`

## Outstanding items before public release

### Complete primary raw data

The supplied GitHub ZIP did not contain the 60 BMG MARS workbooks required by `01_preprocess_plate_data.R`. They must either:

1. be added to `data/raw/`, or
2. be deposited in an external data repository and linked/accessioned in the root README.

Without those files (or a public processed dataset sufficient to start at step 2), the primary manuscript results cannot be reproduced end-to-end from the current repository alone.

### License

No software/data license was selected automatically. Choose and add an appropriate license before making the repository public.

### Runtime validation

The scripts supplied by the user had previously been run during the manuscript analysis. This cleanup environment does not contain R, so the refactored public-path versions could not be executed end-to-end here. The modifications to scripts 01–06 are path/name changes only; scripts 07–09 were refactored to remove the obsolete v6.4 dependency. A final fresh-session run from the repository root is recommended before release.
