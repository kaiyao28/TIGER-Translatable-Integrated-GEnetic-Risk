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

Expected result: the scripts print probabilities, create five figures under
`examples/figures/`, and the tests finish with `All tests passed.`

## 5. Understand the key terms

| Term | Meaning in TIGER |
|---|---|
| PRS | Polygenic risk score |
| RV | Rare variant; separate from APOE |
| `K` | Population disorder prevalence |
| `SP` | Sample prevalence, or case proportion, in the target context |
| `SP_RV` | Background level for the RV update; default 0.50 for the balanced sample calculation |
| `r2_liability` | Variance explained by the PRS on the liability scale |
| PAIR (summary) | Default example conversion using theoretical moments from `K` and liability-scale R² |
| APOE status | Optional genotype column with one value per individual |
| RV status | Optional delimited RV-ID column or mapped logical/0–1 columns |

## 6. Replace the synthetic files carefully

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
