# Disorder-specific and custom RV references

TIGER keeps supplied SCZ and AD reference inputs separate:

```text
inst/extdata/reference/
├── SCZ/
│   ├── SCZ_PTV.xlsx
│   ├── SCZ_MPC2.xlsx
│   └── SCZ_CNV.xlsx
└── AD/
    ├── AD_LOF.xlsx
    ├── AD_REVEL25.xlsx
    ├── AD_REVEL50.xlsx
    ├── AD_REVEL75.xlsx
    └── AD_SRV.xlsx
```

These are literature-derived reference inputs, not universal defaults. Users must confirm
provenance, phenotype, ancestry, frequency definition, effect estimate and reuse
rights before application.

## Load the supplied SCZ references

```r
library(readxl)
library(TIGER)

scz_ptv <- read_tiger_reference(
  tiger_reference_file("SCZ", "SCZ_PTV.xlsx"), disorder = "SCZ", class = "PTV"
)
scz_mpc2 <- read_tiger_reference(
  tiger_reference_file("SCZ", "SCZ_MPC2.xlsx"), disorder = "SCZ", class = "MPC2"
)
scz_cnv <- read_tiger_reference(
  tiger_reference_file("SCZ", "SCZ_CNV.xlsx"), disorder = "SCZ", class = "CNV"
)
```

PTV uses `CAP_freq`, `COP_freq` and `CP_OR`; MPC2 uses `CAM2_freq`,
`COM2_freq` and `CM2_OR`. For CNVs, TIGER calculates the odds ratio from
`Case_freq` and `Control_freq`.

## Load the supplied AD references

```r
ad_lof <- read_tiger_reference(
  tiger_reference_file("AD", "AD_LOF.xlsx"), disorder = "AD", class = "LOF"
)
ad_revel50 <- read_tiger_reference(
  tiger_reference_file("AD", "AD_REVEL50.xlsx"), disorder = "AD", class = "REVEL50"
)
ad_srv <- read_tiger_reference(
  tiger_reference_file("AD", "AD_SRV.xlsx"), disorder = "AD", class = "SRV"
)
```

AD LOF/REVEL tables use gene-level `Case_OR`; the SRV table uses variant-level
`rsID` and `OR`.

The REVEL25, REVEL50 and REVEL75 workbooks represent alternative, potentially
nested thresholds. Select and justify the intended threshold. Do not combine
all threshold files as independent RV classes, because that can count the same
underlying evidence more than once.

## Create a custom reference

Users may supply any disorder-specific evidence after harmonising it to:

| Column | Required meaning |
|---|---|
| `ID` | Unique stable variant, gene or syndrome identifier |
| `Symbol` | Human-readable label |
| `Class` | Variant/evidence class |
| `Case_freq` | Comparable carrier/genotype frequency among cases, 0–1 |
| `Control_freq` | Comparable carrier/genotype frequency among controls, 0–1 |
| `OR` | Positive odds ratio; values below 1 are protective |

```r
library(TIGER)

custom_raw <- read.csv("my_disorder_reference.csv")
custom_reference <- prepare_rv_reference(data.frame(
  ID = custom_raw$variant_id,
  Symbol = custom_raw$gene,
  Class = custom_raw$annotation,
  Case_freq = custom_raw$case_carrier_frequency,
  Control_freq = custom_raw$control_carrier_frequency,
  OR = custom_raw$odds_ratio
), source = "MY_DISORDER")
```

## Checks before use

- Confirm that case and control frequencies use the same unit: carrier,
  genotype or allele frequency.
- Use comparable phenotype definitions, ancestry and ascertainment.
- Check counts, confidence intervals, zero cells and effect uncertainty.
- Do not convert relative risk to an odds ratio by relabelling it.
- Check overlap among genes, variants, CNVs and annotation thresholds.
- Do not combine correlated or nested effects as independent evidence without a
  justified joint model.
- Record the source and transformation applied to every retained row.
- Do not publish identifiable individual-level carrier data.

See [`DATA.md`](DATA.md) for the individual carrier-status schema and
[`REFERENCES.md`](REFERENCES.md) for study attribution.
