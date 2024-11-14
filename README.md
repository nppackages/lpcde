  [![license: GPLv2](https://img.shields.io/badge/license-GPLv2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![CRAN version](https://img.shields.io/cran/v/lpcde?color=7733BB&label=CRAN)](https://cran.r-project.org/web/packages/lpcde/index.html)
[![codecov](https://codecov.io/gh/nppackages/lpcde/branch/main/graph/badge.svg?token=DE4RI272QB)](https://codecov.io/gh/nppackages/lpcde)

# lpcde

The `lpcde` package provides R implementation of bandwidth selection, point estimation and inference procedures for local polynomial conditional distribution and density methods.

This work was supported by the National Science Foundation through grant [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805) and [SES-1947662](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947662), and by the National Institutes of Health through grant [R01 GM072611-16](https://reporter.nih.gov/project-details/10093056).

## Website

https://nppackages.github.io/lpcde

## Queries and Requests

Please email: [lpcde@googlegroups.com](mailto:lpcde@googlegroups.com)

## R Implementation

### Install from source

Install the development version of the package by running

```
devtools::install_github('nppackages/lpcde/R/lpcde')
```

### Install from CRAN
To install/update in R type:

```
install.packages('lpcde')
```

## Contributing
If you have any contributions to help improve or increase the functionality of this package, then please feel free to contribute, by opening a PR or an issue if you have any suggestions.

## Major Updates
Details regarding package updates can be found in [NEWS.md](https://github.com/nppackages/lpcde/blob/main/R/lpcde/NEWS.md).

## Additional Information

- Help: [R Manual](https://cran.r-project.org/web/packages/lpcde/lpcde.pdf), [CRAN repository](https://cran.r-project.org/web/packages/lpcde/index.html).

- A simple replication file: [R-script](R/lpcde_illustration.R).

- Additional simulations using this package: [Technical Paper replication script](https://github.com/nppackages-replication/CCJM_2024_Bernoulli)


## References

### Software and Implementation

- Cattaneo, Chandak, Jansson and Ma (2024): [lpcde: Estimation and Inference for Local Polynomial Conditional Density Estimators](https://nppackages.github.io/references/Cattaneo-Chandak-Jansson-Ma_2024_lpcde.pdf).<br>
Instructional guide.


### Technical and Methodological

- Cattaneo, Chandak, Jansson and Ma (2024): [Boundary Adaptive Local Polynomial Conditional Density Estimators](https://nppackages.github.io/references/Cattaneo-Chandak-Jansson-Ma_2024_Bernoulli.pdf).<br>
Bernoulli 30(4): 3193-3223.<br>
[Supplemental appendix](https://nppackages.github.io/references/Cattaneo-Chandak-Jansson-Ma_2024_Bernoulli--Supplemental.pdf).


<br>

## License
&copy; 2024 lpcde Authors

The contents of this repository are distributed under the MIT license. See below
for details:

```
Copyright (c) 2024 lpcde Authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```