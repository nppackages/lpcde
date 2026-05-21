# Local Polynomial Conditional Density Methods

The `lpcde` package implements bandwidth selection, point estimation, and inference procedures for local polynomial conditional distribution and density methods.

- `lpcde`: local polynomial conditional CDF, PDF, and derivative estimation with pointwise and uniform inference.
- `lpbwcde`: rule-of-thumb bandwidth selection for local polynomial conditional density estimation.

## R Implementation

To install/update in R type:
```
install.packages('lpcde')
```

- Help: [CRAN repository](https://CRAN.R-project.org/package=lpcde).

- Replication: [R illustration](replication/lpcde_illustration.R), [Python illustration](replication/lpcde_illustration.py), [Stata illustration](replication/lpcde_illustration.do), [software article replication](replication/lpcde_replication.R), [comparison illustration](replication/lpcde_comparison.R).

## Python Implementation

To install/update locally from this repository:
```
pip install lpcde
```

- Help: [Python README](Python/lpcde/README.md), [illustration script](replication/lpcde_illustration.py).

- Development usage:
```
from lpcde import lpcde, lpbwcde
```

The Python implementation mirrors the R `lpcde` and `lpbwcde` numerical paths and includes cross-language regression tests against R-generated fixtures.

## Stata Implementation

To install/update in Stata type:
```
net install lpcde, from(https://raw.githubusercontent.com/nppackages/lpcde/main/stata) replace
```

- Commands: `lpcde` for estimation and inference, `lpbwcde` for bandwidth selection.

- Help: [lpcde](stata/lpcde.pdf), [lpbwcde](stata/lpbwcde.pdf).

- Replication: [do-file](replication/lpcde_illustration.do).

The Stata implementation is standalone Stata/Mata code, with numerical smoke tests intended to keep it aligned with the R and Python baselines.


## References

### Software and Implementation

- Cattaneo, Chandak, Jansson and Ma (2025): [lpcde: Estimation and Inference for Local Polynomial Conditional Density Estimators](https://mdcattaneo.github.io/papers/Cattaneo-Chandak-Jansson-Ma_2025_JOSS.pdf).<br>
_Journal of Open Source Software_ 10(107): 7241.<br>
[Companion arXiv version](https://mdcattaneo.github.io/papers/Cattaneo-Chandak-Jansson-Ma_2025_JOSS--arXiv.pdf)


### Technical and Methodological

- Cattaneo, Chandak, Jansson and Ma (2024): [Boundary Adaptive Local Polynomial Conditional Density Estimators](https://nppackages.github.io/references/Cattaneo-Chandak-Jansson-Ma_2024_Bernoulli.pdf).<br>
Bernoulli 30(4): 3193-3223.<br>
[Supplemental appendix](https://nppackages.github.io/references/Cattaneo-Chandak-Jansson-Ma_2024_Bernoulli--Supplemental.pdf).


## Funding

This work was supported in part by the National Science Foundation through grants [SES-1947805](https://www.nsf.gov/awardsearch/show-award/?AWD_ID=1947805), SES-1947662, [DMS-2210561](https://www.nsf.gov/awardsearch/show-award/?AWD_ID=2210561), and [SES-2241575](https://www.nsf.gov/awardsearch/show-award/?AWD_ID=2241575), and by the National Institutes of Health through grant R01 GM072611-16.


<br><br>
