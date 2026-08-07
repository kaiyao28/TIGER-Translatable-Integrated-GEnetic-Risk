# Using TIGER

This guide shows the complete simulation-free application workflow. See
[`DATA.md`](DATA.md) for schemas and [`METHODS.md`](METHODS.md) for method details.

## 1. Prepare the liability-scale PRS

TIGER examples require a correctly centred liability-scale PRS and matching
liability-scale PRS R². Target and population-reference scores must use the same
variants, alleles, effects and centring.

```r
source("R/tiger.R")
load_tiger()

target <- read.csv(tiger_example_file("example_target_prs_observed.csv"))
reference <- read.csv(tiger_example_file("example_reference_prs_observed.csv"))

K <- 0.01
SP_score <- 0.50

converted <- prepare_liability_prs_inputs(
  target_prs_observed = target$PRS_observed,
  reference_prs_observed = reference$PRS_observed,
  K = K,
  SP = SP_score,
  center_on_reference = FALSE
)

target$PRS_liability <- converted$target_prs_liability
r2_liability <- converted$r2_liability
```

Use `center_on_reference = FALSE` when both scores were already centred using
the same population-reference allele frequencies. See
[`LIABILITY_PRS_GUIDE.md`](LIABILITY_PRS_GUIDE.md) for the complete preparation and QC
workflow. A runnable example is in `examples/liability_conversion_example.R`.

## 2. Use one merged individual table

The simplest interface accepts one row per person. The table may contain any
additional columns, including a study group. Only the configured column names
are interpreted by TIGER.

```r
individuals <- data.frame(
  participant_id = c("P001", "P002", "P003"),
  Group = c("Control", "Case", "Case"),
  PRS_liability = c(-0.42, 0.31, 0.18),
  RV_status = c("", "RISK_A", "RISK_A;PROTECT_C"),
  APOE = c("e3/e3", "e3/e4", "e4/e4")
)

rv_reference <- prepare_rv_reference(read.csv(
  tiger_example_file("example_rv_reference.csv")
))
apoe_reference <- prepare_high_impact_reference(read.csv(
  tiger_example_file("example_apoe_reference.csv")
))

results <- tiger_probabilities(
  individuals, K = K, SP = 0.50,
  method = "PAIR (summary)",
  id_col = "participant_id", prs_col = "PRS_liability",
  group_col = "Group", r2_liability = r2_liability,
  include_rv = TRUE, rv_reference = rv_reference,
  rv_status_col = "RV_status", rv_prevalence = 0.50,
  include_apoe = TRUE, apoe_col = "APOE",
  apoe_reference = apoe_reference
)
```

`RV_status` lists reference IDs separated by semicolons, commas, or vertical
bars. Alternatively, use separate 0/1 columns and map them explicitly:

```r
results <- tiger_probabilities(
  individuals_wide, K = K, SP = 0.50,
  r2_liability = r2_liability,
  include_rv = TRUE, rv_reference = rv_reference,
  rv_columns = c(
    RISK_A = "carries_risk_a",
    PROTECT_C = "carries_protect_c"
  )
)
```

The names on the left must match `ID` in the RV reference. The values on the
right name columns in the individual table. The original columns are preserved
and output rows remain aligned to their input IDs.

| Optional layers | Returned probability columns |
|---|---|
| neither | `Probability_PRS` |
| RV | `Probability_PRS`, `Probability_PRS_RV` |
| Generic high-impact component | `Probability_PRS`, `Probability_PRS_HIGH_IMPACT` |
| RV and generic high-impact component | all four conditions, including `Probability_PRS_HIGH_IMPACT_RV` |
| APOE | `Probability_PRS`, `Probability_PRS_APOE` |
| RV and APOE | all four conditions, including `Probability_PRS_APOE_RV` |

All optional layers default to off. `rv_reference` is mandatory when
`include_rv = TRUE`. Use `include_high_impact` for a generic high-impact
component and `include_apoe` for the six-genotype APOE model. Supply the
corresponding genotype column and reference. Beforehand, exclude the separately
modelled variant and a justified LD-aware region from PRS construction and R²
estimation.

## 3. Lower-level, separate-record workflow

The following sections expose the underlying calculations for users who need
to inspect or customise each stage.

### Read the separate input records

```r
individuals <- read.csv(tiger_example_file("example_individuals.csv"))
apoe_status <- read.csv(tiger_example_file("example_apoe_status.csv"))
rv_status <- read.csv(tiger_example_file("example_rv_carriers.csv"))

rv_reference <- prepare_rv_reference(read.csv(
  tiger_example_file("example_rv_reference.csv")
))
apoe_reference <- prepare_high_impact_reference(read.csv(
  tiger_example_file("example_apoe_reference.csv")
))
```

The PRS table contains the individual ID, liability-scale PRS, and optional
plotting group. APOE and RV status remain separate individual-level records.
Reference files contain effect/frequency evidence, not individual status.

### Match APOE status by ID

```r
if (anyDuplicated(individuals$ID) || anyDuplicated(apoe_status$ID)) {
  stop("individual and APOE-status IDs must be unique")
}
apoe_match <- match(individuals$ID, apoe_status$ID)
if (anyNA(apoe_match)) {
  stop("APOE status is missing for one or more individuals")
}
individuals$APOE <- apoe_status$APOE[apoe_match]
```

Never join individual files by assumed row order.

### Calculate the four probability conditions

PAIR (summary) is the worked-example default:

```r
K <- 0.01
SP <- 0.50
r2l <- 0.10
SP_RV <- 0.50

# PRS only
p_prs <- pair_probability_summary(
  individuals$PRS_liability,
  K = K,
  SP = SP,
  r2_liability = r2l
)

# PRS + APOE
p_prs_apoe <- high_impact_method_probability(
  individuals$PRS_liability,
  individuals$APOE,
  apoe_reference,
  K = K, SP = SP, method = "PAIR (summary)",
  r2_liability = r2l
)

# Prepare the separate RV carrier layer
carrier_matrix <- prepare_rv_carrier_matrix(
  rv_status,
  individual_ids = individuals$ID,
  variant_ids = rv_reference$ID
)
rv_effects <- setNames(rv_reference$OR, rv_reference$ID)

# PRS + RV
rv_result <- apply_rv_carriers(
  p_prs,
  carrier_matrix,
  rv_effects,
  prevalence = SP_RV
)

# PRS + APOE + RV
combined_result <- apply_rv_carriers(
  p_prs_apoe,
  carrier_matrix,
  rv_effects,
  prevalence = SP_RV
)
```

`SP_RV = 0.50` defines TIGER's balanced 50:50 case/control sample calculation.
The resulting RV component
is therefore a sample-level overall probability. For a population-level
application using population-level frequency/effect evidence, replace `SP_RV`
with the justified population-level background probability. The chosen level
must match the reference evidence and intended interpretation and must be
reported.

The four results are:

| Condition | Object |
|---|---|
| PRS | `p_prs` |
| PRS + RV | `rv_result$probability_after` |
| PRS + APOE | `p_prs_apoe` |
| PRS + APOE + RV | `combined_result$probability_after` |

APOE is a common high-impact component, not an RV. Either optional component
may be used alone or both may be included.

### Assemble and export a results table

```r
results <- data.frame(
  ID = individuals$ID,
  PRS_liability = individuals$PRS_liability,
  APOE = individuals$APOE,
  RV_count = rv_result$RV_count,
  Probability_PRS = p_prs,
  Probability_PRS_RV = rv_result$probability_after,
  Probability_PRS_APOE = p_prs_apoe,
  Probability_PRS_APOE_RV = combined_result$probability_after
)

write.csv(results, "tiger_results.csv", row.names = FALSE)
```

Inspect row counts, missing values, probability ranges and ID alignment before
using this table. Do not publish individual-level genetic-risk results or
identifiers without the required governance and disclosure controls.

## 4. Plot RV carriers

Unconnected hollow circles show the PRS-only distributions for the two groups.
Larger filled points show RV-adjusted carriers in the matching group colour.
Shape distinguishes one RV from two or more. Colour is not used to encode
damaging or protective direction because that direction is visible from the
probability shift.

![PRS curve with adjusted points only for RV carriers](../examples/figures/rv_carrier_points_v4.png)

```r
individual_group <- individuals$Group
carrier_index <- rv_result$RV_count > 0
carrier_group <- individual_group[carrier_index]
rv_labels_by_id <- tapply(
  rv_status$Variant_ID, rv_status$ID,
  function(x) paste(unique(x), collapse = "; ")
)
rv_labels <- unname(rv_labels_by_id[individuals$ID[carrier_index]])
group_colours <- c("Control" = "#3182BD", "Case" = "#E31A1C")
plot_tiger_rv_carrier_points(
  prs_curve = individuals$PRS_liability,
  probability_curve = p_prs,
  prs_group = individual_group,
  carrier_prs = individuals$PRS_liability[carrier_index],
  carrier_probability_before = p_prs[carrier_index],
  carrier_probability_after = rv_result$probability_after[carrier_index],
  rv_count = rv_result$RV_count[carrier_index],
  rv_labels = rv_labels,
  carrier_group = carrier_group,
  group_colours = group_colours
)
```

Supply one `prs_group` per PRS value and the corresponding `carrier_group` for
each carrier. The same named `group_colours` are applied to hollow PRS
observations and filled adjusted points. RV names are shown when `rv_labels` is supplied. See
[`PLOTTING.md`](PLOTTING.md).

## 5. Plot APOE and combined APOE + RV

APOE-only genotype curves:

![APOE genotype probability curves](../examples/figures/apoe_genotype_curves_v2.png)

```r
prs_sequence <- seq(-4, 4, by = 0.05)
p_sequence <- pair_probability_summary(
  prs_sequence, K = K, SP = SP, r2_liability = r2l
)

plot_tiger_apoe_curves(
  prs_sequence,
  p_sequence,
  apoe_reference,
  K = K, SP = SP, r2_liability = r2l
)
```

Group-specific APOE distributions with final probabilities plotted for RV
carriers:

![APOE curves with adjusted RV-carrier points](../examples/figures/apoe_rv_carrier_points_v1.png)

```r
plot_tiger_apoe_rv_carrier_points(
  prs_curve = individuals$PRS_liability,
  probability_prs = p_prs,
  apoe_reference = apoe_reference,
  carrier_prs = individuals$PRS_liability[carrier_index],
  carrier_apoe = individuals$APOE[carrier_index],
  carrier_probability_before_rv = p_prs_apoe[carrier_index],
  carrier_probability_after_rv =
    combined_result$probability_after[carrier_index],
  rv_count = rv_result$RV_count[carrier_index],
  prs_group = individual_group,
  prs_apoe = individuals$APOE,
  rv_labels = rv_labels,
  carrier_group = carrier_group,
  group_colours = group_colours,
  K = K, SP = SP, r2_liability = r2l
)
```

Observed APOE status is supplied through `prs_apoe`. Each person therefore
contributes one point calculated with the corresponding genotype-specific
`K_g` and `SP_g`. All six genotype distributions are shown together and labelled
directly. Blue and red indicate controls and cases. Filled RV-carrier points
retain the group colour, and shape identifies one or multiple RVs.
All plotting helpers use `tiger_plot_theme()` and shared APOE/RV conventions.

## Run the complete example

```bash
Rscript examples/example_data_and_plots.R
Rscript tests/run_tests.R
```
