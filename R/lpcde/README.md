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

- Replication: [R-script](https://github.com/nppackages/lpcde/blob/main/replication/lpcde_illustration.R), [software article replication](https://github.com/nppackages/lpcde/blob/main/replication/lpcde_replication.R), [comparison illustration](https://github.com/nppackages/lpcde/blob/main/replication/lpcde_comparison.R), [Python illustration](https://github.com/nppackages/lpcde/blob/main/replication/lpcde_illustration.py), [Stata illustration](https://github.com/nppackages/lpcde/blob/main/replication/lpcde_illustration.do).

- Development version:
```
devtools::install_github('nppackages/lpcde/R/lpcde')
```

Basic usage:
```
model1 <- lpcde::lpcde(x_data = x_data, y_data = y_data,
                       y_grid = y_grid, x = 0, bw = 1)

model2 <- lpcde::lpbwcde(y_data = y_data, x_data = x_data,
                         x = 0, y_grid = y_grid)
```

Standard R methods including `coef`, `confint`, `plot`, `print`, `summary`, and `vcov` are available for fitted `lpcde` objects. The `summary`, `print`, and `coef` methods are available for `lpbwcde` bandwidth-selection objects.

## References

### Software and Implementation

- Cattaneo, Chandak, Jansson and Ma (2025): [lpcde: Estimation and Inference for Local Polynomial Conditional Density Estimators](https://mdcattaneo.github.io/papers/Cattaneo-Chandak-Jansson-Ma_2025_JOSS.pdf).<br>
_Journal of Open Source Software_ 10(107): 7241.<br>
[Companion arXiv version](https://mdcattaneo.github.io/papers/Cattaneo-Chandak-Jansson-Ma_2025_JOSS--arXiv.pdf).

### Technical and Methodological

- Cattaneo, Chandak, Jansson and Ma (2024): [Boundary Adaptive Local Polynomial Conditional Density Estimators](https://nppackages.github.io/references/Cattaneo-Chandak-Jansson-Ma_2024_Bernoulli.pdf).<br>
Bernoulli 30(4): 3193-3223.<br>
[Supplemental appendix](https://nppackages.github.io/references/Cattaneo-Chandak-Jansson-Ma_2024_Bernoulli--Supplemental.pdf).

## Funding

This work was supported in part by the National Science Foundation through grants [SES-1947805](https://www.nsf.gov/awardsearch/show-award/?AWD_ID=1947805), SES-1947662, [DMS-2210561](https://www.nsf.gov/awardsearch/show-award/?AWD_ID=2210561), and [SES-2241575](https://www.nsf.gov/awardsearch/show-award/?AWD_ID=2241575), and by the National Institutes of Health through grant R01 GM072611-16.
