# Method references

Please cite every method used in an analysis.

## BPC

Uffelmann E, Major Depressive Disorder Working Group of the Psychiatric
Genomics Consortium, Schizophrenia Working Group of the Psychiatric Genomics
Consortium, et al. Estimating disorder probability based on polygenic
prediction using the BPC approach. *Nature Communications*. 2025;16:8443.
https://doi.org/10.1038/s41467-025-62929-x

Authors' implementation (BPC v1.0.1):
https://doi.org/10.5281/zenodo.15721085

## Pain absolute-scale conversion

Pain O, Gillett AC, Austin JC, Folkersen L, Lewis CM. A tool for translating
polygenic scores onto the absolute scale using summary statistics. *European
Journal of Human Genetics*. 2022;30:339–348.
https://doi.org/10.1038/s41431-021-01028-z

## Common- and rare-variant probability framework

Escott-Price V, Schmidt KM. Probability of Alzheimer's disease based on common
and rare genetic variants. *Alzheimer's Research & Therapy*. 2021;13:140.
https://doi.org/10.1186/s13195-021-00884-7

## Adaptations in this repository

- `pain_probability(corrected = TRUE)` exposes a population-prevalence moment
  adaptation. `corrected = FALSE` gives the original sample-prevalence version.
- `pair_probability_summary()` supplies theoretical liability-scale
  case/control moments to the population-logistic framework and is labelled
  **PAIR (summary)**.
- `pair_probability_sample()` accepts case/control PRS moments estimated from a
  labelled sample and is labelled **PAIR (sample)**. In-sample moments give a
  descriptive comparator; prospective use requires an independent calibration
  dataset.
- `apply_rv_probability()` supports damaging and protective RV effects.

Distinguish adaptations from original methods in publications and cite the
source equations together with the software version used.
