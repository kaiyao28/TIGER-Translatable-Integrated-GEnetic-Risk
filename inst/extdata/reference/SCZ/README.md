# SCZ reference inputs

This folder contains the SCZ reference inputs distributed with TIGER:

- `SCZ_PTV.xlsx`: gene PTV carrier frequencies and `CP_OR`;
- `SCZ_MPC2.xlsx`: gene MPC > 2 carrier frequencies and `CM2_OR`; and
- `SCZ_CNV.xlsx`: CNV/syndrome case and control frequencies.

Load them with `read_tiger_reference()` from `R/reference_data.R`. The CNV odds
ratio is derived from case/control frequencies rather than copied from the
workbook's relative-risk column.

These files are literature-derived references, not universal SCZ defaults. Check
the source publication, phenotype, ancestry, counts, uncertainty and redistribution terms before
reuse. Do not assume that all retained rows or classes are mutually independent.
See `../../../../docs/REFERENCE_DATA.md`.
