# Input-data schemas

Users are required to prepare their own RV reference. See
[`REFERENCE_DATA.md`](REFERENCE_DATA.md) for the required schema.

## Harmonised schema

| Field | Type | Definition |
|---|---|---|
| `RV_IDs` | character | Required carrier-matching key; one value per reference row |
| `OR` | numeric | Positive odds ratio; required unless calculated from frequencies |
| `Symbol` | character | Optional human-readable label; defaults to `RV_IDs` |
| `Class` | character | Optional PTV, MPC > 2, CNV or other category |
| `Case_freq` | numeric | Case carrier frequency used when calculating a missing OR |
| `Control_freq` | numeric | Control carrier frequency used when calculating a missing OR |
| `Population_freq` | numeric | Optional replacement for a missing control frequency |

Column names are matched case-insensitively. Each row must contain a positive
OR or enough frequency information to calculate one. When `Control_freq` is
missing, `Population_freq` is used as its proxy. The harmonised output records
the provenance in `OR_source` and `Control_freq_source`.
For imported legacy tables, `ID`, `RV_ID`, and `Variant_ID` are accepted aliases
for `RV_IDs` and are normalised internally.

For a custom table, TIGER calculates a missing odds ratio from comparable case
and control carrier frequencies:

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
| `High_impact_genotype` or `APOE_genotype` | Genotype label matching the corresponding individual column |
| `Case_freq` | Genotype frequency among cases |
| `Control_freq` | Genotype frequency among controls |

`Genotype` is accepted as an imported-table alias for either genotype field and
is normalised internally.

Case and control frequencies must each sum to one. For APOE, observed genotype
frequencies are preferred. `apoe_genotype_reference()` can derive the six
genotype frequencies from named e2/e3/e4 case and control allele frequencies,
but this makes a Hardy-Weinberg equilibrium assumption. See
[`AD_APOE_GUIDE.md`](AD_APOE_GUIDE.md) for application and PRS-overlap guidance.
For one biallelic common high-impact variant, use observed frequencies for the
0-, 1-, and 2-effect-allele genotypes when available.
`biallelic_genotype_reference()` can derive them from separate case and control
effect-allele frequencies under HWE. See
[`HIGH_IMPACT_VARIANTS.md`](HIGH_IMPACT_VARIANTS.md).

## Lower-level separate APOE-status schema

APOE status may be maintained in a separate source table before constructing
the merged input:

| Field | Definition |
|---|---|
| `ID` | Unique individual identifier matching the prepared PRS input |
| `APOE_genotype` | One genotype label such as `e3/e3`, `e3/e4` or `e4/e4` |

Join APOE status to the PRS table by `ID`, never by row position. Check for
duplicated IDs, unmatched individuals, unexpected genotype labels and missing
values before calling `high_impact_method_probability()`. The genotype labels
must match the `APOE_genotype` values in the APOE case/control reference.

## Merged individual-table schema

`tiger_probabilities()` can read PRS, optional group information, APOE status,
and RV status from one table. One row must represent one unique individual.

| Example field | Definition |
|---|---|
| `ID` | Unique individual identifier |
| `PRS_liability` | Finite liability-scale PRS |
| `Group` (optional) | Group retained unchanged for analysis or plotting |
| `High_impact_genotype` (optional) | Value matching `High_impact_genotype` in the generic reference, such as 0/1/2 effect alleles |
| `APOE_genotype` (optional) | Value matching `APOE_genotype` in the six-genotype APOE reference |
| `RV_IDs` (optional) | Delimited values matching `RV_IDs` in the RV reference |

`High_impact_genotype` and `APOE_genotype` illustrate alternative uses of one high-impact layer;
do not add both in the same calculation. Calculation-input names can be
configured through `id_col`, `prs_col`, `apoe_col` or `high_impact_col`, and
`rv_status_col`. The optional `Group` column is preserved unchanged for later
analysis or plotting. Instead of `RV_IDs`, `rv_columns` can map
RV-reference IDs to separate logical or 0/1 columns. An empty status means no
carried RV. Unknown RV IDs and unknown APOE genotypes produce errors.

## Lower-level separate carrier-status schema

For direct use of `prepare_rv_carrier_matrix()` and `apply_rv_carriers()`, RV
status may instead be kept separately from the RV effect reference. This
lower-level interface accepts a long table with:

| Field | Definition |
|---|---|
| `ID` | Individual identifier matching the PRS input |
| `Variant_ID` | Identifier matching `RV_IDs` in the RV reference |
| `Carrier` | `TRUE`/`FALSE` or `1`/`0` |

Use `prepare_rv_carrier_matrix()` to convert this table to the matrix required
by `apply_rv_carriers()`. A full table may contain every individual/variant
pair. A presence-only table may contain carrier rows only, but then pass the
complete `individual_ids` and `variant_ids` arguments explicitly so omitted
pairs are intentionally interpreted as non-carriers.

The files under `inst/extdata/example/` are synthetic demonstrations only:

- `example_individuals.csv`: IDs, liability-scale PRS, and optional plotting
  groups for separate group-specific PRS distributions. The synthetic example
  contains 200 observations for each APOE genotype, balanced as 100 controls
  and 100 cases;
- `example_apoe_status.csv`: a separate individual-level ID/APOE-genotype
  record joined to the PRS input by ID;
- `example_target_prs_observed.csv`: synthetic target scores before liability
  conversion;
- `example_reference_prs_observed.csv`: synthetic matching population-reference
  scores used for conversion and liability-scale R² estimation;
- `example_rv_reference.csv`: three hypothetical RV effects;
- `example_rv_carriers.csv`: presence-only carrier rows, including people with
  one RV and more than one RV;
- `example_apoe_reference.csv`: hypothetical case/control genotype frequencies.

Copyable starting schemas are provided under `inst/extdata/templates/` for the
individual table, RV reference, APOE reference, and generic high-impact
reference. Locate them after installation with `tiger_template_file()`.

Do not place identifiable individual-level data in a public repository.

## PRS input prerequisite

The `PRS_liability` field in the example individual table is already on the
liability scale. Users must prepare their own target and population-reference
PRSs consistently and perform the required liability conversion before using
TIGER. See [`LIABILITY_PRS_GUIDE.md`](LIABILITY_PRS_GUIDE.md) and the runnable
`examples/liability_conversion_example.R`.

## Data boundary and provenance

Do not add identifiable individual data or third-party reference datasets to
the package. Keep empirical inputs outside the repository and document their
citations, provenance, definitions, and reuse terms in the applied analysis.
