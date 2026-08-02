# AD reference inputs

This folder contains the AD reference inputs distributed with TIGER:

- `AD_LOF.xlsx`: gene-level loss-of-function evidence;
- `AD_REVEL25.xlsx`, `AD_REVEL50.xlsx`, `AD_REVEL75.xlsx`: alternative
  REVEL-threshold gene sets; and
- `AD_SRV.xlsx`: individually identified rare variants.

Load them with `read_tiger_reference()` from `R/reference_data.R`.

The three REVEL thresholds may contain nested/overlapping evidence. Select and
justify the intended threshold; do not stack all three as independent effects.
These files are literature-derived references, not universal AD defaults. Confirm source publication,
phenotype, ancestry, frequency definition, uncertainty and redistribution terms
before reuse. See `../../../../docs/REFERENCE_DATA.md`.
