# Attribution and modification notice

TIGER is under development in preparation for publication. Pull requests are
welcome, but every change must be reviewed and approved by the TIGER
maintainers before merging. This notice records method attribution, not a
running history of intermediate development edits.

TIGER is a reorganized and extended implementation prepared in 2026 for
applying polygenic-score probability conversions and optional independent
rare-variant and separately modelled common high-impact probability layers.

The BPC equations and liability-scale input workflow are based on the BPC method
and public guidance by Uffelmann et al.:

- Publication: https://doi.org/10.1038/s41467-025-62929-x
- Authors' software archive: https://doi.org/10.5281/zenodo.15721085
- Authors' repository: https://github.com/euffelmann/bpc

The Pain quantile conversion is based on Pain et al.:

- https://doi.org/10.1038/s41431-021-01028-z

The population-logistic and independent damaging-RV framework is based on
Escott-Price and Schmidt:

- https://doi.org/10.1186/s13195-021-00884-7

The implementation reorganizes the cited methods behind a consistent interface
and adds documented TIGER integrations for protective RVs, common high-impact
genotypes, APOE and simulation-independent visualization. Methodological
adaptations are described in the function documentation and REFERENCES.md.

This notice provides attribution and identifies modifications. It does not
alter the terms of the GNU General Public License supplied in `LICENSE`.
