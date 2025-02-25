# LPCDE
The `lpcde` package provides R implementation for bandwidth selection,
point estimation and inference for local polynomial conditional distribution
and density methods.

This work was supported by the the National Science Foundation through grants
SES-1947805, SES-1947662, DMS-2210561, and SES-2241575,
and from the National Institute of Health (R01 GM072611-16).

## Website
https://nppackages.github.io/lpcde

## R Implementation
To install/update in R type:
```
install.packages('lpcde')
```
Here we outline some of the simple use cases of the package. For a given dataset, 
comprising of a response variable (`y_data`) and a collection of regressors (`x_data`),
the `lpcde` function estimates the conditional distribution (or derivatives thereof)
at some prescribe value of `x`. By default a locally quadratic polynomial in y 
and locally linear polynomial in x and Epanechnikov kernel are used on  some prescribed 
values of y, `y_grid`, with a bandwidth `bw`.
This procedure is executed with the following command:
```
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1)
```
Standard R methods, `coef`, `confint`, `vcov`, `print`, `plot` and `summary`, 
can be used on objects returned from the `lpcde` function to understand the output.
For example, the summary command applied to the object generated from this function,
`summary(model1)`, will output a table that provides the point estimates of the 
conditional density, number of data points used in estimation, standard error, 
and a robust 95% confidence interval. The confidence level, robustness of point 
and standard error estimates, estimates of uniform confidence bands can all be 
altered by providing the correct inputs to the optional variables for the `lpcde`
function.

By default, the function provides estimates according to the original
formulation of the estimator (see references below for all technical details).
If a constrained density estimate that is non-negative and integrates to one is 
desired (see references below), 
the flags `nonneg` and `normalize` can be turned on:
```
model_reg = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1,
                        nonneg=TRUE, normalize=TRUE)
```

The package also provides additional an function for automating the bandwidth 
selection, `lpbwcde`. By default this function computes the rule-of-thumb MSE 
optimal bandwidth for the conditional PDF with locally quadratic polynomial in y 
and locally linear polynomial in x and Epanechnikov kernel on the implied support
of Y.
```
model2 = lpbwcde(y_data = y_data, x_data = x_data, x = 0, y_grid = y_grid)

```
The estimated bandwidth from this function can be used as bandwidth input to
`lpcde` directly by using the option `bwselect` in `lpcde` to specify bandwidth
selection type instead of running `lpbwcde` first.
R methods like `summary`, `print` and `coef` can be applied to the objects 
returned by `lpbwcde` for understanding the output.

Additional replication scripts showcasing the full functionality of the package 
along with references for a full description of the methods along with
theoretical guarantees can be found on the GitHub repository `lpcde`.

## References
### Software and Implementation
- Cattaneo, Chandak, Jansson and Ma (2025):
lpcde: Estimation and Inference for Local Polynomial Conditional Density Estimators
<br>
_JOSS_, forthcoming.

### Technical and Methodological
-  Cattaneo, Chandak, Jansson and Ma (2024): 
Boundary Adaptive Local Polynomial Conditional Density Estimators
<br>
Bernoulli 30(4): 3193-3223.
