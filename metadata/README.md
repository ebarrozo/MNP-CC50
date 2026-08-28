# Plate-layout metadata

`key_plate_layout.xlsx` records the experimental 96-well plate map.

Rows A–H:

- A: PS Small
- B: PE Small
- C: PVC Small
- D: PS Large
- E: PE Large
- F: PVC Large
- G: Mix Small
- H: Mix Large

Experimental wells are columns 1–9 in three technical-replicate dose blocks. Controls are:

- A10–A12: PBS
- B10–B12: Methanol
- C10–C12: Media

## Dose correction

The historical plate key contains working-suspension labels of 1, 10, and 100 µg/mL. Fifty microliters of working suspension was added to 200 µL already present in each experimental well, giving final well concentrations of **0.2, 2, and 20 µg/mL**. The analysis scripts use the corrected final concentrations and do not derive dose values from the historical labels in the key.
