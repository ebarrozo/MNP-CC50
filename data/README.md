# Data

## `raw/`

Place the 60 original BMG MARS EthD-1 plate-reader `.xlsx` exports here for end-to-end reproduction of the primary analysis.

Expected design:

- 5 cell lines
- 3 independent experiments
- 4 exposure durations
- 60 workbooks total

Expected filename pattern:

```text
EthD-1 MARS_<time>hr[_btmoptic]_<cell line>_<replicate>.xlsx
```

Examples:

```text
EthD-1 MARS_6hr_btmoptic_Jurkat_3.xlsx
EthD-1 MARS_24hr_THP-1_2.xlsx
```

The public preprocessing script expects the original BMG MARS worksheet named `End point` and the plate matrix at the format used by the source exports.

## `no_cell/`

Contains the MNP-only, no-cell EthD-1 optical-interference experiment used for Supplementary Table S2.

The no-cell experiment was performed without the wash steps used immediately before EthD-1 acquisition in cell-containing wells. It is therefore used to characterize potential particle-associated fluorescence rather than as a direct subtraction control.
