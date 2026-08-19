# Preparing an RV reference

TIGER does not distribute literature-derived, study-specific, or individual-
level genetic data. Files under `inst/extdata/example/` are synthetic or
hypothetical and exist only to demonstrate the interface.

Users must select, document, and supply an appropriate RV effect reference for
their disorder, population, phenotype definition, and intended probability
scale.

## Recommended schema

Use one row per independently modelled RV or RV burden:

| Column | Meaning |
|---|---|
| `RV_IDs` | Unique value matching the individual `RV_IDs` field |
| `OR` | Positive odds ratio, unless calculated from frequencies |
| `Symbol` | Optional human-readable annotation |
| `Class` | Optional variant or burden class |
| `Case_freq` | Optional case carrier frequency |
| `Control_freq` | Optional control carrier frequency |
| `Population_freq` | Optional proxy when control frequency is unavailable |

```r
rv_reference <- prepare_rv_reference(data.frame(
  RV_IDs = c("RISK_A", "PROTECT_C"),
  Symbol = c("GENE_A", "GENE_C"),
  Class = c("Hypothetical_risk", "Hypothetical_protective"),
  OR = c(8.0, 0.4)
), source = "HYPOTHETICAL_EXAMPLE")
```

The identifiers and values above are hypothetical. Replace them with reviewed
evidence and record its source outside the package.

## Calculating a missing odds ratio

If `OR` is missing, TIGER calculates it from comparable carrier frequencies:

```text
OR = [Case_freq / (1 - Case_freq)] /
     [Control_freq / (1 - Control_freq)]
```

When `Control_freq` is missing, `Population_freq` may be used as a proxy. The
prepared table records these choices in `OR_source` and
`Control_freq_source`. This substitution must be scientifically justified.

## Importing an Excel table

`read_tiger_reference()` optionally maps user-supplied Excel files that follow
the supported SCZ or AD source layouts:

```r
reference <- read_tiger_reference(
  path = "path/to/my_reference.xlsx",
  disorder = "AD",
  class = "LOF"
)
```

For any other layout, import the file normally, rename or map its columns, and
pass the resulting data frame to `prepare_rv_reference()`.

## Checks before use

- Confirm comparable phenotype, ancestry, and ascertainment definitions.
- Confirm whether frequencies describe alleles, genotypes, or carriers.
- Review confidence intervals, zero cells, small counts, and effect uncertainty.
- Avoid double-counting overlapping variants, genes, CNVs, or nested classes.
- Ensure separately modelled high-impact regions are excluded from the PRS.
- Preserve citations and provenance in the user's analysis records.
