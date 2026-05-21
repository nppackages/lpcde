# LPCDE

The `lpcde` Python package implements local polynomial conditional distribution and density estimation with bandwidth selection.

- `lpcde`: local polynomial conditional CDF, PDF, and derivative estimation with pointwise and uniform inference quantities.
- `lpbwcde`: rule-of-thumb bandwidth selection for local polynomial conditional density estimation.

Website: [https://nppackages.github.io/lpcde/](https://nppackages.github.io/lpcde/).

Source code: [https://github.com/nppackages/lpcde](https://github.com/nppackages/lpcde).

## Authors

Matias D. Cattaneo (<matias.d.cattaneo@gmail.com>)

Rajita Chandak (<rajita.chandak@gmail.com>)

Michael Jansson (<michael.jansson.berkeley@gmail.com>)

Xinwei Ma (<xinweima.pku@gmail.com>)

## Installation

To install locally from the repository:
```
pip install -e Python/lpcde
```

## Usage

```
import numpy as np
from lpcde import confint, lpbwcde, lpcde, plot, summary, vcov

x_data = np.random.normal(size=500)
y_data = np.random.normal(loc=x_data, scale=1)
y_grid = np.linspace(-1, 1, 5)

fit = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=0.5)
bw = lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid=y_grid)
summary(fit)
confint(fit)
```

Replication examples are available in `replication/lpcde_illustration.py` and `replication/lpcde_replication.py`.

## Dependencies

- numpy
- pandas
- scipy
- matplotlib

## References

### Software and Implementation

- Cattaneo, Chandak, Jansson and Ma (2025): [lpcde: Estimation and Inference for Local Polynomial Conditional Density Estimators](https://mdcattaneo.github.io/papers/Cattaneo-Chandak-Jansson-Ma_2025_JOSS.pdf).<br>
_Journal of Open Source Software_ 10(107): 7241.

### Technical and Methodological

- Cattaneo, Chandak, Jansson and Ma (2024): [Boundary Adaptive Local Polynomial Conditional Density Estimators](https://nppackages.github.io/references/Cattaneo-Chandak-Jansson-Ma_2024_Bernoulli.pdf).<br>
Bernoulli 30(4): 3193-3223.
