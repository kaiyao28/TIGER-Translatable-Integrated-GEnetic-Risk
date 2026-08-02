# TIGER: Translatable Integrated Genetic Risk framework

TIGER is a simulation-free R framework for estimating absolute binary-disorder
probabilities from polygenic risk, rare variants and separately modelled common
high-impact variants such as APOE.

```text
                         ┌─ PRS only
PRS ─> absolute          ├─ PRS + RV
      probability ───────┼─ PRS + APOE / common high-impact variant
                         └─ PRS + APOE / high-impact variant + RV
```

APOE and RV are separate optional components. APOE is not treated as an RV.

## Start here

1. Prepare a correctly centred liability-scale PRS and matching liability-scale
   PRS R².
2. Supply separate individual records for PRS, APOE status and RV carrier
   status.
3. Supply justified case/control evidence for APOE and RV effects.
4. Choose and justify population prevalence `K`, sample prevalence `SP`, and the RV
   calculation level (`SP_RV`; default 0.50 for a balanced sample).
5. Calculate PRS only, PRS + RV, PRS + APOE, or the combined probability.

The worked-example defaults are:

```r
K <- 0.01
SP <- 0.50
r2_liability <- 0.10
SP_RV <- 0.50
```

PAIR (summary) is the default example method. These values are instructional,
not universal disease parameters.

## Quick example

New to TIGER? Start with [`GETTING_STARTED.md`](docs/GETTING_STARTED.md). Core
calculations use R 4.1+; figures require `ggplot2` and supplied Excel references
require `readxl`.

If using RStudio, open `TIGER.Rproj` first. This sets the repository as the
working project. Then run the complete synthetic workflow in
`examples/example_data_and_plots.R` from top to bottom. The shorter code below
is organised in cumulative sections and assumes earlier sections have run.

Install the development version directly from GitHub:

```r
install.packages("remotes") # once, if needed
remotes::install_github("OWNER/TIGER", build_vignettes = FALSE)
library(TIGER)
```

Replace `OWNER` with the repository owner after publication. Package functions
remain plain, reviewable R source under `R/`; installation does not hide the
implementation.

From the `TIGER` directory:

```bash
Rscript examples/liability_conversion_example.R
Rscript examples/example_data_and_plots.R
Rscript tests/run_tests.R
```

The examples use synthetic inputs under `inst/extdata/example/` and produce all four
probability conditions plus the documented plots. No files outside this
repository are required.

## Three levels of implementation

The levels are cumulative. Begin by converting PRS to probabilities, then add
RV information, and finally add APOE or another separately modelled common
high-impact variant. Complete data-loading and ID-matching code is provided in
[`USAGE.md`](docs/USAGE.md).

### Level 1: PRS-to-probability conversion

TIGER includes the three conversion approaches evaluated by the framework:
BPC, GenoPred, and PAIR. PAIR is available with theoretical summary-derived
moments or with externally estimated case/control sample moments.

```r
# BPC
p_bpc <- bpc_probability(
  individuals$PRS_liability, K = K, SP = SP,
  r2_liability = r2_liability
)

# GenoPred
r2_observed <- liability_to_observed_r2(r2_liability, K = K, SP = SP)
p_genopred <- genopred_probability(
  individuals$PRS_liability,
  reference_prs_liability,
  r2_observed = r2_observed,
  K = K,
  SP = SP
)

# PAIR (summary), the worked-example default
p_prs <- pair_probability_summary(
  individuals$PRS_liability, K, SP, r2_liability
)
```

`PAIR (sample)` is available through `pair_probability_sample()` when
case/control PRS moments have been estimated in an appropriate independent
calibration sample.

### Level 2: PRS + rare variants

Apply the RV layer to the selected PRS probability. The PRS curve is common to
all externally defined groups. Group differences are shown by the adjusted
carrier points rather than by separate PRS curves.

Example inputs:

**Individual PRS input**

| ID | PRS_liability |
| --- | ---: |
| P001 | -3.50 |
| P002 | -3.00 |

**RV effect reference**

| ID | Symbol | Class | Case_freq | Control_freq | OR |
| --- | --- | --- | ---: | ---: | ---: |
| RISK_A | GENE_A | PTV | 0.010 | 0.002 | 8.0 |
| PROTECT_C | GENE_C | Protective | 0.003 | 0.008 | 0.4 |

**Individual RV-carrier input**

| ID | Variant_ID | Carrier |
| --- | --- | ---: |
| P005 | RISK_A | 1 |
| P005 | RISK_B | 1 |

Repeated IDs indicate multiple RVs. Omitted IDs are non-carriers in this
presence-only format. See the complete synthetic
[PRS](inst/extdata/example/example_individuals.csv),
[RV reference](inst/extdata/example/example_rv_reference.csv), and
[RV-carrier](inst/extdata/example/example_rv_carriers.csv) files and the
[schema](docs/DATA.md#individual-carrier-status-schema).

```r
rv_result <- apply_rv_carriers(
  p_prs, carrier_matrix, odds_ratios, prevalence = SP_RV
)

# TRUE for each person carrying at least one RV; FALSE for non-carriers.
carrier_index <- rv_result$RV_count > 0

# Make an explicit carrier-only table for plotting.
carriers <- individuals[carrier_index, , drop = FALSE]
carriers$Probability_before_RV <- p_prs[carrier_index]
carriers$Probability_after_RV <-
  rv_result$probability_after[carrier_index]
carriers$RV_count <- rv_result$RV_count[carrier_index]
rv_labels_by_id <- tapply(
  rv_carriers$Variant_ID, rv_carriers$ID,
  function(x) paste(unique(x), collapse = "; ")
)
carriers$RV_label <- unname(rv_labels_by_id[carriers$ID])
carriers$Carrier_group <- carrier_groups$Carrier_group[
  match(carriers$ID, carrier_groups$ID)
]
group_colours <- c("Group A" = "#3182BD", "Group B" = "#E31A1C")

plot_tiger_rv_carrier_points(
  prs_curve = prs_sequence,
  probability_curve = p_sequence,
  carrier_prs = carriers$PRS_liability,
  carrier_probability_before = carriers$Probability_before_RV,
  carrier_probability_after = carriers$Probability_after_RV,
  rv_count = carriers$RV_count,
  rv_labels = carriers$RV_label,
  carrier_group = carriers$Carrier_group,
  group_colours = group_colours
)
```

`carrier_index` is `TRUE` for carriers. The resulting table keeps all plotting
inputs in the same individual order. `drop = FALSE` preserves a data frame when
only one carrier is present.

The curve is PRS only. Points are shown only for RV carriers. A circle denotes
one RV and a triangle denotes two or more RVs. The example uses two external
groups to demonstrate point colours. Point colour never represents damaging or
protective RV direction.

See the [figure-customisation guide](docs/PLOTTING.md#rv-carrier-points) to
show RV names, separate nearby labels with guide lines, or change point size,
opacity, shape, border, group colours and legends.

![PRS curve with adjusted points only for RV carriers](examples/figures/rv_carrier_points_v4.png)

### Level 3: PRS + APOE or another high-impact variant + RV

First recalculate the selected PRS conversion using genotype-specific
population and sample prevalence. Then apply the same RV layer to the
APOE-updated probabilities.

Example inputs, joined to the PRS table by `ID`:

**Individual APOE-status input**

| ID | APOE |
| --- | --- |
| P001 | e2/e3 |
| P003 | e3/e4 |

**APOE case/control reference**

| Genotype | Case_freq | Control_freq |
| --- | ---: | ---: |
| e2/e2 | 0.0016 | 0.0064 |
| e4/e4 | 0.0900 | 0.0169 |

The complete reference must contain all six genotypes, with case and control
frequencies each summing to one. See the synthetic
[APOE-status](inst/extdata/example/example_apoe_status.csv) and
[APOE-reference](inst/extdata/example/example_apoe_reference.csv) files and the
[schema](docs/DATA.md#individual-apoe-status-schema).

```r
p_prs_apoe <- high_impact_method_probability(
  individuals$PRS_liability, individuals$APOE, apoe_reference,
  K = K, SP = SP, method = "PAIR (summary)", r2_liability = r2_liability
)

plot_tiger_apoe_curves(
  prs_sequence, p_sequence, apoe_reference,
  K = K, SP = SP, r2_liability = r2_liability
)

combined_result <- apply_rv_carriers(
  p_prs_apoe, carrier_matrix, odds_ratios, prevalence = SP_RV
)

# Add the combined-model probabilities to the same carrier-only table.
carriers$Probability_APOE_before_RV <- p_prs_apoe[carrier_index]
carriers$Probability_APOE_after_RV <-
  combined_result$probability_after[carrier_index]

plot_tiger_apoe_rv_carrier_points(
  prs_curve = prs_sequence,
  probability_prs = p_sequence,
  apoe_reference = apoe_reference,
  carrier_prs = carriers$PRS_liability,
  carrier_apoe = carriers$APOE,
  carrier_probability_before_rv = carriers$Probability_APOE_before_RV,
  carrier_probability_after_rv = carriers$Probability_APOE_after_RV,
  rv_count = carriers$RV_count,
  rv_labels = carriers$RV_label,
  carrier_group = carriers$Carrier_group,
  group_colours = group_colours,
  K = K, SP = SP, r2_liability = r2_liability
)
```

The APOE-only figure shows one updated PRS curve per genotype. In the combined
figure, those curves stay common across external groups. Coloured points show
the final group-specific RV-carrier observations, while shape distinguishes one
from multiple RVs. APOE is modelled separately and is not treated as an RV.
The same [RV point and label controls](docs/PLOTTING.md#rv-carrier-points) apply
to this combined figure.

![APOE genotype probability curves](examples/figures/apoe_genotype_curves_v2.png)

![APOE curves with adjusted points only for RV carriers](examples/figures/apoe_rv_carrier_points_v1.png)

The complete runnable code—including creation of `prs_sequence` and carrier
subsets—is in [`USAGE.md`](docs/USAGE.md) and
`examples/example_data_and_plots.R`.

## Documentation

| Guide                                            | Contents                                                 |
| ------------------------------------------------ | -------------------------------------------------------- |
| [`GETTING_STARTED.md`](docs/GETTING_STARTED.md) | Installation, first run, glossary and common errors      |
| [`USAGE.md`](docs/USAGE.md)                     | Complete application code and plotting examples          |
| [`METHODS.md`](docs/METHODS.md)                 | PAIR summary/sample, BPC, GenoPred, APOE and RV methods  |
| [`PLOTTING.md`](docs/PLOTTING.md)               | Point styling, group colours, legends and ggplot editing |
| [`DATA.md`](docs/DATA.md)                       | PRS, APOE, RV-reference and individual-status schemas    |
| [`REFERENCE_DATA.md`](docs/REFERENCE_DATA.md)   | Separate SCZ/AD references and custom-reference creation |
| [`BPC_INPUT_GUIDE.md`](docs/BPC_INPUT_GUIDE.md) | Liability-scale PRS preparation and quality control      |
| [`AD_APOE_GUIDE.md`](docs/AD_APOE_GUIDE.md)     | APOE-region exclusion and AD/APOE application            |
| [`REFERENCES.md`](docs/REFERENCES.md)           | Method publications and adaptations                      |

## Repository structure

```text
R/                 loader, probability, reference, high-impact, RV and plotting functions
inst/extdata/example/      synthetic input-format examples
inst/extdata/reference/    separate SCZ and AD reference inputs
examples/          runnable application examples and rendered figures
tests/             dependency-light checks
```

## Scope

TIGER does not run GWAS, calculate PRSs from genotype files, select scientific
inputs, establish clinical utility or replace external validation. Users must
justify ancestry, prevalence, SP, score construction, effect evidence,
overlap and independence assumptions. See [`USAGE.md`](docs/USAGE.md) and the focused guides
before application.

## Development

TIGER is under development in preparation for publication. Cloning, testing,
discussion and reviewed pull requests are welcome. See
[`CONTRIBUTING.md`](docs/CONTRIBUTING.md), [`NOTICE.md`](docs/NOTICE.md),
[`CITATION.cff`](CITATION.cff) and the GPL v3-or-later [`LICENSE`](LICENSE).
