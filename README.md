<p align="center">
  <img src="docs/assets/tiger-readme-header-final.png" alt="TIGER — Translatable Integrated GEnetic Risk" width="100%">
</p>

TIGER is an R package for converting a liability-scale polygenic risk score
(PRS) into a disorder probability and optionally incorporating rare variants
(RVs) and one separately modelled common high-impact genetic component.

```text
Liability-scale PRS ──> estimated disorder probability
                         ├─ PRS only
                         ├─ PRS + RV
                         ├─ PRS + common high-impact component
                         │    ├─ one biallelic variant: 0/1/2 effect alleles
                         │    └─ APOE: six e2/e3/e4 genotypes
                         └─ PRS + common high-impact component + RV
```

The selected method converts PRS to probability. TIGER then optionally adds
RVs and/or one common high-impact component. APOE is a six-genotype example of
the high-impact layer, not an RV.

`K` is population prevalence and `SP` is the target-sample case proportion.
`rv_prevalence` is the probability prior for the RV update, not RV frequency.
See [METHODS.md](docs/METHODS.md) for definitions and scale selection.

## Installation

```r
install.packages("remotes") # once, if needed
remotes::install_github(
  "kaiyao28/TIGER-Translatable-Integrated-Genetic-Risk",
  build_vignettes = FALSE
)
library(TIGER)
```

Core calculations require R 4.1 or later. Optional plotting functions require
`ggplot2`, with `ggrepel` used when available to improve RV-label placement.
The `read_tiger_reference()` helper requires `readxl` when importing Excel
reference files.

## Probability methods

The four methods below convert **PRS to probability**. TIGER applies the
optional high-impact and RV updates if required.

- `"BPC"`
- `"GenoPred"`
- `"PAIR (summary)"`—theoretical moments from liability-scale R²
- `"PAIR (sample)"`—case/control PRS moments from a suitable sample

The examples use `"PAIR (summary)"`. See [METHODS.md](docs/METHODS.md) for
required inputs, assumptions, and formulas.

## 1. Prepare the liability-scale PRS

All TIGER probability methods require a correctly centred liability-scale PRS.
Target and population-reference scores must use the same variants, effect
alleles, effect estimates, and reference-frequency centring.

```r
converted <- prepare_liability_prs_inputs(
  target_prs_observed = target$PRS_observed,
  reference_prs_observed = reference$PRS_observed,
  K = 0.01,
  SP = 0.50,
  center_on_reference = FALSE
)

individuals$PRS_liability <- converted$target_prs_liability
reference_prs_liability <- converted$reference_prs_liability
r2_liability <- converted$r2_liability
```

Use `center_on_reference = FALSE` when both input scores were already centred
using the same population-reference allele frequencies. Do not independently
standardise the target and reference PRSs.

If reference-derived liability-scale R² is unavailable, a compatible
independently evaluated or leave-one-out observed-scale R² may be converted:

```r
r2_liability <- observed_to_liability_r2(
  r2_observed = reported_r2_observed,
  K = 0.01,
  SP = SP_evaluation
)
```

`SP_evaluation` is the case proportion in the sample where the observed-scale
R² was estimated. See the [liability-scale PRS guide](docs/LIABILITY_PRS_GUIDE.md)
for score construction, centring, R² alternatives, and quality control.

## 2. Prepare one individual table

One row represents one person. Only `ID` and `PRS_liability` are required.
Group, RV, and high-impact genotype columns are optional. `APOE` below is one
example of a high-impact genotype column.

| ID   | PRS_liability | Group (optional) | RV_status (optional) | High_impact (optional) | APOE (optional) |
| ---- | ------------: | ---------------- | -------------------- | ---------------------- | --------------- |
| P001 |         -0.42 | Control          |                      | 0                      | e3/e3           |
| P002 |          0.31 | Case             | RISK_A               | 1                      | e3/e4           |
| P003 |          0.18 | Case             | RISK_A;PROTECT_C     | 2                      | e4/e4           |

Column names can be specified when calling `tiger_probabilities()`. RV status
may be a delimited list of RV-reference IDs or separate logical/0–1 columns.
`Group` is optional plotting information only and does not affect any
probability calculation.
`High_impact` illustrates a single 0/1/2 variant, whereas `APOE` illustrates
the six-genotype APOE model. They are alternative examples; use only the column
for the component being modelled.

When RVs are included, provide one reference row per RV. Its `ID` must match
the value used in `RV_status`:

| ID        | Symbol | Class      | Case_freq | Control_freq |  OR |
| --------- | ------ | ---------- | --------: | -----------: | --: |
| RISK_A    | GENE_A | PTV        |     0.010 |        0.002 | 8.0 |
| PROTECT_C | GENE_C | Protective |     0.003 |        0.008 | 0.4 |

When APOE is included, provide mutually exclusive genotype frequencies among
cases and controls. Each frequency column must sum to one:

| Genotype | Case_freq | Control_freq |
| -------- | --------: | -----------: |
| e2/e2    |    0.0016 |       0.0064 |
| e2/e3    |    0.0528 |       0.1264 |
| e2/e4    |    0.0240 |       0.0208 |
| e3/e3    |    0.4356 |       0.6241 |
| e3/e4    |    0.3960 |       0.2054 |
| e4/e4    |    0.0900 |       0.0169 |

These values are illustrative. Use compatible study references. See
[DATA.md](docs/DATA.md) for schemas and
[REFERENCE_DATA.md](docs/REFERENCE_DATA.md) for preparation and provenance.
The same genotype-reference schema supports a single 0/1/2 high-impact
variant; see [HIGH_IMPACT_VARIANTS.md](docs/HIGH_IMPACT_VARIANTS.md).

## 3. Calculate probabilities

`tiger_probabilities()` is the main application function. The example below
uses APOE as the common high-impact component and then adds RVs.

Before enabling a high-impact component, exclude that variant and a justified
LD-aware region from PRS construction and R² estimation. Apply the same
exclusion to target and reference PRSs and document the genome build and
excluded interval or variant set. TIGER cannot verify this preprocessing step.

```r
results <- tiger_probabilities(
  data = individuals,
  K = 0.01,
  SP = 0.50,
  method = "PAIR (summary)",
  id_col = "ID",
  prs_col = "PRS_liability",
  group_col = "Group",
  r2_liability = r2_liability,
  include_rv = TRUE,
  rv_reference = rv_reference,
  rv_status_col = "RV_status",
  rv_prevalence = 0.50,
  include_apoe = TRUE,
  apoe_col = "APOE",
  apoe_reference = apoe_reference
)
```

Here `rv_prevalence = 0.50` gives a balanced-sample RV update. Population-level
use requires a compatible reference and justified population prior.

The input columns and row order are preserved. TIGER always returns
`Probability_PRS` and adds clearly labelled columns for enabled components:

| Enabled component                               | Additional probability columns                                                   |
| ----------------------------------------------- | -------------------------------------------------------------------------------- |
| RV                                              | `Probability_PRS_RV`                                                           |
| Generic or three-genotype high-impact component | `Probability_PRS_HIGH_IMPACT`, and with RV: `Probability_PRS_HIGH_IMPACT_RV` |
| Six-genotype APOE component                     | `Probability_PRS_APOE`, and with RV: `Probability_PRS_APOE_RV`               |

RV count and damaging/protective component columns are also returned when RVs
are enabled.

### PRS only

```r
results_prs <- tiger_probabilities(
  individuals,
  K = 0.01, SP = 0.50,
  method = "PAIR (summary)",
  r2_liability = r2_liability
)
```

`results_prs$Probability_PRS` contains the PRS-derived probability. The optional
`plot_tiger_prs_methods()` example below compares the available methods.

<img src="examples/figures/prs_method_comparison_v1.png" alt="PRS probability-method comparison" width="650">

### PRS + RV

```r
results_rv <- tiger_probabilities(
  individuals,
  K = 0.01, SP = 0.50,
  method = "PAIR (summary)",
  r2_liability = r2_liability,
  include_rv = TRUE,
  rv_reference = rv_reference,
  rv_status_col = "RV_status",
  rv_prevalence = 0.50
)
```

`results_rv` contains `Probability_PRS` and `Probability_PRS_RV`, allowing the
probability before and after RV incorporation to be compared directly. The
optional example uses `plot_tiger_rv_carrier_points()`.

<img src="examples/figures/rv_carrier_points_v4.png" alt="PRS and RV-adjusted carrier probabilities" width="650">

### PRS + one common high-impact variant

Represent one biallelic variant by mutually exclusive genotypes with 0, 1, or
2 effect alleles. Observed case/control genotype frequencies are preferred.
When only effect-allele frequencies are available, the helper below assumes
HWE separately in cases and controls:

```r
# Replace these illustrative frequencies with values for your variant.
single_variant_reference <- biallelic_genotype_reference(
  case_effect_allele_frequency = 0.20,
  control_effect_allele_frequency = 0.10
)

results_single_variant <- tiger_probabilities(
  individuals,
  K = 0.01, SP = 0.50,
  method = "PAIR (summary)",
  r2_liability = r2_liability,
  include_rv = TRUE,
  rv_reference = rv_reference,
  rv_status_col = "RV_status",
  rv_prevalence = 0.50,
  include_high_impact = TRUE,
  high_impact_col = "Variant_genotype",
  high_impact_reference = single_variant_reference
)
```

`Variant_genotype` must contain `0`, `1`, or `2` effect alleles. Replace the
illustrative frequencies with study values. The helper assumes HWE within
cases and controls; observed genotype frequencies are preferred.

<img src="examples/figures/single_high_impact_variant_rv_carrier_points_v1.png" alt="Single common high-impact variant distributions and RV-adjusted carriers" width="650">

Full assumptions and reference options are in
[HIGH_IMPACT_VARIANTS.md](docs/HIGH_IMPACT_VARIANTS.md).

### PRS + APOE

```r
# Start with TIGER's illustrative six-genotype reference.
apoe_reference <- tiger_default_apoe_reference()

# Inspect it, then replace the frequencies with a suitable reference when available.
apoe_reference

results_apoe <- tiger_probabilities(
  individuals,
  K = 0.01, SP = 0.50,
  method = "PAIR (summary)",
  r2_liability = r2_liability,
  include_apoe = TRUE,
  apoe_col = "APOE",
  apoe_reference = apoe_reference
)
```

The bundled frequencies demonstrate the workflow only. Replace them with an
appropriate reference when available.

<img src="examples/figures/apoe_genotype_curves_v2.png" alt="APOE genotype-specific probability distributions" width="650">

The same interface accepts another mutually exclusive multi-variant genotype
system. See [HIGH_IMPACT_VARIANTS.md](docs/HIGH_IMPACT_VARIANTS.md).

### PRS + APOE + RV

```r
results_combined <- tiger_probabilities(
  individuals,
  K = 0.01, SP = 0.50,
  method = "PAIR (summary)",
  r2_liability = r2_liability,
  include_rv = TRUE,
  rv_reference = rv_reference,
  rv_status_col = "RV_status",
  rv_prevalence = 0.50,
  include_apoe = TRUE,
  apoe_col = "APOE",
  apoe_reference = apoe_reference
)
```

`results_combined` contains all four probability conditions. The optional
`plot_tiger_apoe_rv_carrier_points()` example shows RV-adjusted carrier
probabilities relative to the APOE genotype distributions.

<img src="examples/figures/apoe_rv_carrier_points_v1.png" alt="APOE genotype distributions and RV-adjusted carriers" width="650">

## Plotting is optional

The returned table can be used with any plotting software. TIGER's optional
helpers return ordinary `ggplot` objects. See
[PLOTTING.md](docs/PLOTTING.md) for interpretation and customisation, or run
[example_data_and_plots.R](examples/example_data_and_plots.R).

## Run the examples

From the repository root:

```bash
Rscript examples/liability_conversion_example.R
Rscript examples/example_data_and_plots.R
Rscript tests/run_tests.R
```

The example data are synthetic and require no external individual-level data.

## Documentation

| Guide                                                  | Purpose                                              |
| ------------------------------------------------------ | ---------------------------------------------------- |
| [GETTING_STARTED.md](docs/GETTING_STARTED.md)           | Installation, first run, glossary, and common errors |
| [USAGE.md](docs/USAGE.md)                               | Complete high- and lower-level application workflows |
| [LIABILITY_PRS_GUIDE.md](docs/LIABILITY_PRS_GUIDE.md)   | Liability-scale PRS preparation and R² inputs       |
| [DATA.md](docs/DATA.md)                                 | Individual, RV, and high-impact reference schemas    |
| [METHODS.md](docs/METHODS.md)                           | Probability methods, APOE update, and RV integration |
| [PLOTTING.md](docs/PLOTTING.md)                         | Optional plotting and customisation                  |
| [REFERENCE_DATA.md](docs/REFERENCE_DATA.md)             | Supplied SCZ/AD references and custom references     |
| [AD_APOE_GUIDE.md](docs/AD_APOE_GUIDE.md)               | APOE-region exclusion and AD/APOE application        |
| [HIGH_IMPACT_VARIANTS.md](docs/HIGH_IMPACT_VARIANTS.md) | Single common high-impact variant application        |
| [REFERENCES.md](docs/REFERENCES.md)                     | Method publications and adaptations                  |

## Scope

TIGER does not run GWAS, construct PRSs from genotype files, select scientific
inputs, establish clinical utility, or replace external validation. Users must
justify ancestry, population prevalence `K`, sample prevalence `SP`, PRS and
reference compatibility, variant independence, and overlap between separately
modelled components.

TIGER is under development in preparation for publication. See
[CONTRIBUTING.md](docs/CONTRIBUTING.md), [NOTICE.md](docs/NOTICE.md),
[CITATION.cff](CITATION.cff), and the GPL v3-or-later [LICENSE](LICENSE).
