# Input-data schemas

Supplied disorder references are separated under `inst/extdata/reference/SCZ/` and
`inst/extdata/reference/AD/`. See [`REFERENCE_DATA.md`](REFERENCE_DATA.md) for the file
list, disorder-specific column mappings, loading functions and instructions for
creating a custom reference. These workbooks are reference inputs, not
dependencies or universal defaults.

## Harmonised schema

| Field | Type | Definition |
|---|---|---|
| `ID` | character | Optional stable record identifier |
| `Symbol` | character | Gene, syndrome or variant label |
| `Class` | character | Optional PTV, MPC > 2, CNV or other category |
| `Case_freq` | numeric | Proportion of cases carrying the variant |
| `Control_freq` | numeric | Proportion of controls carrying the variant |
| `OR` | numeric | Positive odds ratio; OR < 1 is protective |

For a custom table, an odds ratio may be calculated from comparable case and
control carrier frequencies:

```r
OR <- (Case_freq / (1 - Case_freq)) /
      (Control_freq / (1 - Control_freq))
```

Example custom harmonisation:

```r
library(TIGER)

raw <- read.csv("my_reference.csv")
reference <- prepare_rv_reference(raw, source = "MY_DISORDER")
```

## Interpretation cautions

- Confirm whether each frequency is an allele, genotype or carrier frequency.
- Confirm comparable ancestry and phenotype definitions in cases and controls.
- Do not interpret an odds ratio as a probability.
- Review zero cells, confidence intervals and small counts.
- Avoid double-counting overlapping CNVs, genes or variant classes.
- Ensure separately modelled high-effect variants are excluded from the PRS
  when required by the scientific model.

## Common high-impact genotype schema

`prepare_high_impact_reference()` accepts mutually exclusive genotype rows:

| Field | Definition |
|---|---|
| `Genotype` | Genotype label matching the target data |
| `Case_freq` | Genotype frequency among cases |
| `Control_freq` | Genotype frequency among controls |

Case and control frequencies must each sum to one. For APOE, observed genotype
frequencies are preferred. `apoe_genotype_reference()` can derive the six
genotype frequencies from named e2/e3/e4 case and control allele frequencies,
but this makes a Hardy-Weinberg equilibrium assumption. See
[`AD_APOE_GUIDE.md`](AD_APOE_GUIDE.md) for application and PRS-overlap guidance.

## Individual APOE-status schema

Keep APOE status in a separate individual-level table:

| Field | Definition |
|---|---|
| `ID` | Unique individual identifier matching the prepared PRS input |
| `APOE` | One genotype label such as `e3/e3`, `e3/e4` or `e4/e4` |

Join APOE status to the PRS table by `ID`, never by row position. Check for
duplicated IDs, unmatched individuals, unexpected genotype labels and missing
values before calling `high_impact_method_probability()`. The genotype labels
must match the `Genotype` values in the APOE case/control reference.

## Individual carrier-status schema

Individual RV status should be supplied separately from the RV effect
reference. TIGER accepts a long table with:

| Field | Definition |
|---|---|
| `ID` | Individual identifier matching the PRS input |
| `Variant_ID` | Identifier matching `ID` in the RV reference |
| `Carrier` | `TRUE`/`FALSE` or `1`/`0` |

Use `prepare_rv_carrier_matrix()` to convert this table to the matrix required
by `apply_rv_carriers()`. A full table may contain every individual/variant
pair. A presence-only table may contain carrier rows only, but then pass the
complete `individual_ids` and `variant_ids` arguments explicitly so omitted
pairs are intentionally interpreted as non-carriers.

The files under `inst/extdata/example/` are synthetic demonstrations only:

- `example_individuals.csv`: IDs and liability-scale PRS;
- `example_apoe_status.csv`: a separate individual-level ID/APOE genotype
  record joined to the PRS input by ID;
- `example_target_prs_observed.csv`: synthetic target scores before liability
  conversion;
- `example_reference_prs_observed.csv`: synthetic matching population-reference
  scores used for conversion and liability-scale R² estimation;
- `example_rv_reference.csv`: three illustrative RV effects;
- `example_rv_carriers.csv`: presence-only carrier rows, including people with
  one RV and more than one RV;
- `example_carrier_groups.csv`: two synthetic external groups used only to
  demonstrate carrier-point colours; and
- `example_apoe_reference.csv`: illustrative case/control genotype frequencies.

Do not place identifiable individual-level data in a public repository.

## PRS input prerequisite

The `PRS_liability` field in the example individual table is already on the
liability scale. Users must prepare their own target and population-reference
PRSs consistently and perform the required liability conversion before using
TIGER. See [`BPC_INPUT_GUIDE.md`](BPC_INPUT_GUIDE.md) and the runnable
`examples/liability_conversion_example.R`.

## Redistribution and provenance

The software licence does not automatically grant rights to third-party data.
Each distributed workbook must carry its original study citation, provenance
and applicable reuse terms. Where redistribution is not permitted, distribute a
download script or a synthetic schema example instead of the source workbook.

Never commit individual-level genotype or phenotype data.
