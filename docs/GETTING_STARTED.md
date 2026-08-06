# Getting started with TIGER

This page is for readers applying TIGER for the first time. The supplied data
are synthetic examples; begin with them before substituting study data.

## 1. Requirements

- R 4.1 or later;
- `ggplot2` for figures; and
- `readxl` for the supplied SCZ/AD Excel references.

Install the optional packages once:

```r
install.packages(c("ggplot2", "readxl"))
```

The probability calculations otherwise use base R and `stats`.

## 2. Open the repository folder

Start R in the folder containing `README.md`, `R/`, `inst/` and `examples/`.
Check the location with:

```r
getwd()
file.exists("R/tiger.R")
```

The second command should return `TRUE`.

## 3. Load TIGER

```r
source("R/tiger.R")
load_tiger()
```

If working from another directory, supply the repository path:

```r
load_tiger(root = "/path/to/TIGER")
```

## 4. Run the supplied examples

From a terminal in the repository folder:

```bash
Rscript examples/liability_conversion_example.R
Rscript examples/example_data_and_plots.R
Rscript tests/run_tests.R
```

Expected result: the scripts print probabilities, create seven figures under
`examples/figures/`, and the tests finish with `All tests passed.`

## 5. Understand the key terms

| Term | Meaning in TIGER |
|---|---|
| PRS | Polygenic risk score |
| RV | Rare variant; separate from APOE |
| `K` | Population disorder prevalence |
| `SP` | Sample prevalence, or case proportion, in the target context |
| `SP_RV` | Disease-probability prior for the RV update; 0.50 gives the balanced-sample calculation and is not RV allele frequency |
| `r2_liability` | Variance explained by the PRS on the liability scale |
| PAIR (summary) | Default example conversion using theoretical moments from `K` and liability-scale R² |
| APOE status | Optional genotype column with one value per individual |
| RV status | Optional delimited RV-ID column or mapped logical/0–1 columns |

## 6. Start with PRS only

The smallest useful call needs an ID, a liability-scale PRS, `K`, `SP`, and the
method-specific PRS input:

```r
individuals_prs <- data.frame(
  ID = c("P001", "P002", "P003"),
  PRS_liability = c(-0.50, 0.00, 0.50)
)

results_prs <- tiger_probabilities(
  individuals_prs,
  K = 0.01,
  SP = 0.50,
  method = "PAIR (summary)",
  r2_liability = 0.10
)

results_prs[c("ID", "Probability_PRS")]
```

Use this first to confirm that the PRS inputs and probability conversion work.
The optional RV and common high-impact layers can then be added independently
or together. APOE is one supported high-impact component.

## 7. Add RV and a high-impact layer

This self-contained example uses APOE as the high-impact component. It defines
the individual data and both required references before calculating all four
probability conditions:

```r
K <- 0.01
SP <- 0.50
r2_liability <- 0.10

individuals <- data.frame(
  ID = sprintf("P%03d", 1:6),
  PRS_liability = c(-0.80, -0.25, 0.05, 0.30, 0.65, 1.00),
  Group = c("Control", "Control", "Control", "Case", "Case", "Case"),
  RV_status = c("", "PROTECT_C", "RISK_A", "", "RISK_B",
                "RISK_A;PROTECT_C"),
  APOE = c("e2/e2", "e2/e3", "e2/e4", "e3/e3", "e3/e4", "e4/e4")
)

rv_reference <- data.frame(
  ID = c("RISK_A", "RISK_B", "PROTECT_C"),
  Symbol = c("GENE_A", "GENE_B", "GENE_C"),
  Class = c("PTV", "Damaging_missense", "Protective"),
  Case_freq = c(0.010, 0.015, 0.003),
  Control_freq = c(0.002, 0.005, 0.008),
  OR = c(8.0, 4.0, 0.4)
)

apoe_reference <- data.frame(
  Genotype = c("e2/e2", "e2/e3", "e2/e4", "e3/e3", "e3/e4", "e4/e4"),
  Case_freq = c(0.0016, 0.0528, 0.0240, 0.4356, 0.3960, 0.0900),
  Control_freq = c(0.0064, 0.1264, 0.0208, 0.6241, 0.2054, 0.0169)
)

# PRS + RV only
results_rv <- tiger_probabilities(
  individuals, K = K, SP = SP,
  method = "PAIR (summary)", r2_liability = r2_liability,
  include_rv = TRUE, rv_reference = rv_reference,
  rv_status_col = "RV_status", rv_prevalence = SP
)

# PRS + APOE only
results_apoe <- tiger_probabilities(
  individuals, K = K, SP = SP,
  method = "PAIR (summary)", r2_liability = r2_liability,
  include_high_impact = TRUE, high_impact_col = "APOE",
  high_impact_reference = apoe_reference
)

# Full PRS + RV + APOE model
results <- tiger_probabilities(
  data = individuals,
  K = K,
  SP = SP,
  method = "PAIR (summary)",
  id_col = "ID",
  prs_col = "PRS_liability",
  group_col = "Group",
  r2_liability = r2_liability,
  include_rv = TRUE,
  rv_reference = rv_reference,
  rv_status_col = "RV_status",
  rv_prevalence = SP,
  include_high_impact = TRUE,
  high_impact_col = "APOE",
  high_impact_reference = apoe_reference
)

results[c(
  "ID", "Probability_PRS", "Probability_PRS_RV",
  "Probability_PRS_HIGH_IMPACT", "Probability_PRS_HIGH_IMPACT_RV"
)]
```

`Group` is retained for plotting but does not affect the calculations. Empty
`RV_status` values mean that the individual carries none of the listed RVs.
Multiple RV IDs are separated by a semicolon. The reference values above are
illustrative and must be replaced with suitable study references.

## 8. Replace the synthetic files carefully

Work through this order:

1. Prepare target and population-reference PRSs on the required scale.
2. Prepare one individual table with unique IDs, liability-scale PRSs, and an
   optional `Group` column for plotting.
3. If using APOE, add one genotype column and prepare a genotype-frequency
   reference.
4. If using RVs, add a delimited RV-status column or mapped logical/0–1 columns
   and prepare a harmonised RV effect reference.
5. Match any source records by ID before constructing the merged table.
6. Replace and justify `K`, `SP`, `SP_RV` and liability-scale R².
7. Call `tiger_probabilities()` and inspect the requested probability columns.

Use [`USAGE.md`](USAGE.md) for the complete code, [`DATA.md`](DATA.md) for file
schemas, and [`REFERENCE_DATA.md`](REFERENCE_DATA.md) for supplied or custom
SCZ/AD references.

## Common errors

### “TIGER files were not found”

R is not running from the repository folder. Check `getwd()` or pass the correct
folder to `load_tiger(root = ...)`.

### “requires the ggplot2/readxl package”

Install the missing package with `install.packages()` and restart the command.

### Missing or duplicated IDs

Do not reorder files manually. Make IDs unique and use the checked ID-matching
code in `USAGE.md`.

### Unknown APOE genotype or RV ID

The individual-status label does not match the corresponding reference table.
Standardise labels before calculation; do not silently drop unmatched records.

### Implausible probabilities

Stop and recheck PRS scaling, `K`, `SP`, `SP_RV`, R², ancestry, frequencies,
effect direction, overlap and independence assumptions. Do not interpret the
result clinically until it has been externally validated.
