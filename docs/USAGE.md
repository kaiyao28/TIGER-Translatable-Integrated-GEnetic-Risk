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
SP <- 0.50

converted <- prepare_bpc_inputs(
  target_prs_observed = target$PRS_observed,
  reference_prs_observed = reference$PRS_observed,
  K = K,
  SP = SP,
  center_on_reference = FALSE
)

target$PRS_liability <- converted$target_prs_liability
r2_liability <- converted$r2_liability
```

Use `center_on_reference = FALSE` when both scores were already centred using
the same population-reference allele frequencies. See
[`BPC_INPUT_GUIDE.md`](BPC_INPUT_GUIDE.md) for the complete preparation and QC
workflow. A runnable example is in `examples/liability_conversion_example.R`.

## 2. Read the separate input records

```r
individuals <- read.csv(tiger_example_file("example_individuals.csv"))
apoe_status <- read.csv(tiger_example_file("example_apoe_status.csv"))
rv_status <- read.csv(tiger_example_file("example_rv_carriers.csv"))
carrier_groups <- read.csv(tiger_example_file("example_carrier_groups.csv"))

rv_reference <- prepare_rv_reference(read.csv(
  tiger_example_file("example_rv_reference.csv")
))
apoe_reference <- prepare_high_impact_reference(read.csv(
  tiger_example_file("example_apoe_reference.csv")
))
```

APOE status and RV status are individual-level records separate from the PRS.
Reference files contain effect/frequency evidence, not individual status.

## 3. Match APOE status by ID

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

## 4. Calculate the four probability conditions

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

## 5. Assemble and export a results table

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

## 6. Plot RV carriers

The line is PRS only. Points are shown only for RV carriers. Shape distinguishes
one RV from two or more. All points use one neutral fill unless an external
`carrier_group` is supplied. Colour is not used to encode damaging or protective
direction because that direction is already visible from the probability shift.

![PRS curve with adjusted points only for RV carriers](../examples/figures/rv_carrier_points_v4.png)

```r
prs_sequence <- seq(-4, 4, by = 0.05)
p_sequence <- pair_probability_summary(
  prs_sequence, K = K, SP = SP, r2_liability = r2l
)
carrier_index <- rv_result$RV_count > 0
carrier_group <- carrier_groups$Carrier_group[
  match(individuals$ID[carrier_index], carrier_groups$ID)
]
rv_labels_by_id <- tapply(
  rv_status$Variant_ID, rv_status$ID,
  function(x) paste(unique(x), collapse = "; ")
)
rv_labels <- unname(rv_labels_by_id[individuals$ID[carrier_index]])
group_colours <- c("Group A" = "#3182BD", "Group B" = "#E31A1C")
plot_tiger_rv_carrier_points(
  prs_curve = prs_sequence,
  probability_curve = p_sequence,
  carrier_prs = individuals$PRS_liability[carrier_index],
  carrier_probability_before = p_prs[carrier_index],
  carrier_probability_after = rv_result$probability_after[carrier_index],
  rv_count = rv_result$RV_count[carrier_index],
  rv_labels = rv_labels,
  carrier_group = carrier_group,
  group_colours = group_colours
)
```

To compare externally defined groups, supply one `carrier_group` per carrier
and optional named `group_colours`. RV names are shown when `rv_labels` is
supplied. Point opacity defaults to 1 and can be reduced with
`rv_point_alpha`. See [`PLOTTING.md`](PLOTTING.md).

## 7. Plot APOE and combined APOE + RV

APOE-only genotype curves:

![APOE genotype probability curves](../examples/figures/apoe_genotype_curves_v2.png)

```r
plot_tiger_apoe_curves(
  prs_sequence,
  p_sequence,
  apoe_reference,
  K = K, SP = SP, r2_liability = r2l
)
```

APOE curves with final probabilities plotted only for RV carriers:

![APOE curves with adjusted RV-carrier points](../examples/figures/apoe_rv_carrier_points_v1.png)

```r
plot_tiger_apoe_rv_carrier_points(
  prs_curve = prs_sequence,
  probability_prs = p_sequence,
  apoe_reference = apoe_reference,
  carrier_prs = individuals$PRS_liability[carrier_index],
  carrier_apoe = individuals$APOE[carrier_index],
  carrier_probability_before_rv = p_prs_apoe[carrier_index],
  carrier_probability_after_rv =
    combined_result$probability_after[carrier_index],
  rv_count = rv_result$RV_count[carrier_index],
  rv_labels = rv_labels,
  carrier_group = carrier_group,
  group_colours = group_colours,
  K = K, SP = SP, r2_liability = r2l
)
```

All plotting helpers use `tiger_plot_theme()` and shared APOE/RV visual
conventions.

## Run the complete example

```bash
Rscript examples/example_data_and_plots.R
Rscript tests/run_tests.R
```
