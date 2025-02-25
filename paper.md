---
title: 'lpcde: Estimation and Inference for Local Polynomial Conditional Density Estimators'
tags:
- R
- statistics
- density estimation
- kernel methods
- local polynomials
date: "21 August 2024"
output:
  pdf_document: default
  html_document:
    df_print: paged
authors:
- name: Matias D. Cattaneo
  orcid: 0000-0003-0493-7506
  affiliation: 1
- name: Rajita Chandak
  orcid: 0009-0006-4289-2520
  corresponding: yes
  affiliation: 2
- name: Michael Jansson
  orcid: 0000-0003-4678-7518
  affiliation: 3
- name: Xinwei Ma
  orcid: 0000-0001-8827-9146
  affiliation: 4
bibliography: paper.bib
link-citations: yes
affiliations:
- name: Department of Operations Research and Financial Engineering, Princeton University,
    USA
  index: 1
- name: Institute of Mathematics, EPFL, Switzerland
  index: 2
- name: Department of Economics, University of California, Berkeley, USA
  index: 3
- name: Department of Economics, University of California, San Diego, USA
  index: 4
---

# Summary

Conditional cumulative distribution functions (CDFs), conditional probability
density functions (PDFs), and derivatives thereof, are important parameters of
interest in statistics, econometrics, and other data science disciplines. The
package `lpcde` implements new estimation and inference methods for conditional
CDFs, conditional PDFs, and derivatives thereof, employing the kernel-based
local polynomial smoothing approach introduced in @CCJM_2024_Bernoulli.

The package `lpcde` offers data-driven (pointwise and uniform) estimation and
inference methods for conditional CDFs, conditional PDFs, and derivatives
thereof, which are automatically valid at both interior and boundary points of
the support of the outcome and conditioning variables. For point estimation, the
package offers mean squared error optimal bandwidth selection and associated
optimal mean square and uniform point estimators. For inference, the package
offers valid confidence intervals and confidence bands based on robust
bias-correction techniques [@Calonico-Cattaneo-Farrell_2018_JASA;
@Calonico-Cattaneo-Farrell_2022_Bernoulli]. Finally, these statistical
procedures can be easily used for visualization and graphical presentation of
smooth estimates of conditional CDFs, conditional PDFs, and derivative thereof,
with custom `ggplot` [@ggplot2] commands built for the package.

This package is currently the only open source implementation of an estimator
offering boundary adaptive, data-driven conditional density estimation with
robust bias-corrected pointwise confidence interval and uniform confidence band
constructions, providing users with statistical tools to better understand the
reliability of their empirical analysis. A detailed tutorial, replication files,
and other information on how to use the package can be found in the [GitHub
repository](https://github.com/nppackages/lpcde) and through the [CRAN
repository](https://cran.r-project.org/web/packages/lpcde/index.html). See also
the `lpcde` package website (https://nppackages.github.io/lpcde/) and the
companion arXiv article [@CCJM_2024_lpcde] for additional methodological
information and numerical results.

# Statement of need

@Wand-Jones_1995_Book, @Fan-Gijbels_1996_Book, @simonoff2012smoothing, and
@scott2015multivariate give textbook introductions to kernel-based density and
local polynomial estimation and inference methods. The core idea underlying the
estimator implemented in `lpcde` is to use kernel-based local polynomial
smoothing methods to construct an automatic boundary adaptive estimator for
conditional CDFs, conditional PDFs, and derivatives thereof. The estimator
implemented in this package consists of two steps. The first step estimates the
conditional distribution function using standard local polynomial regression
methods, and the second step applies local polynomial smoothing to the
(non-smooth) local polynomial conditional CDF estimate from the first step to
obtain a smooth estimate of the conditional CDF, conditional PDF, and
derivatives thereof.

A distinct advantage of this estimation method over existing ones is its
boundary adaptivity for a possibly unknown compact support of the data.
Furthermore, the estimator has a simple closed form representation, which leads
to easy and fast implementation. Unlike other boundary adaptive procedures, the
estimation procedures implemented in the package `lpcde` do not require
pre-processing of data, and thus avoid the challenges of hyper-parameter tuning:
only one bandwidth parameter needs to be selected for implementation. See
@CCJM_2024_Bernoulli and @CCJM_2024_lpcde for more details.

# Comparing and contrasting existing toolsets

The package `lpcde` contributes to a small set of open source statistical
software packages implementing estimation and inference methods for conditional
CDF, conditional PDF, and derivatives thereof. More specifically, we identified
three `R` packages, `hdrcde` [@hdrcde], `haldensify` [@haldensify], and `np`
[@np], and one `Python` package, `cde` [@rothfuss2019conditional], which provide
related methodology. There are no open source `Stata` packages that implement
comparable estimation and inference methods. The table below summarizes some of
the main differences between those other packages and `lpcde`. Notably, `lpcde`
is the only package available that provides both pointwise and uniform
uncertainty quantification, in addition to producing boundary adaptive mean
square and uniformly optimal point estimates via data-driven, optimal tuning
parameter selection. Furthermore, the `lpcde` package produces proper
conditional density estimates that are non-negative and integrate to one. These
features are unique contributions of the package to the `R` toolkit and, more
broadly, to the open source statistical community.

| Package  | Programming language | CDF/Derivative estimation | Regularized density | Valid at boundary | Standard error | Valid inference | Confidence bands | Bandwidth selection |
|--------|:------:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|
| `hdrcde` | R      | x | x | x | x | x | x | $\checkmark$ |
| `np`     | R      | x | x | x | $\checkmark$ | x | x | $\checkmark$ |
| `haldensify` | R      | x | x | x | $\checkmark$ | x | x | $\checkmark$ |
| `cde`    | Python | x | x | x | x | x  | x | $\checkmark$ |
| `lpcde`  | R      | $\checkmark$ | $\checkmark$ | $\checkmark$ | $\checkmark$ | $\checkmark$  | $\checkmark$  | $\checkmark$   |

# Acknowledgements

Cattaneo gratefully acknowledges financial support from the National Science Foundation through grants SES-1947805, DMS-2210561, and SES-2241575, and from the National Institute of Health (R01 GM072611-16).
Jansson gratefully acknowledges financial support from the National Science Foundation through grant SES-1947662.

# References
