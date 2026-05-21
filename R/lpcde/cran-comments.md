## Release summary

This is a major release updating the package to version 1.0.0.

This release:

* updates the maintainer to Matias D. Cattaneo;
* updates package metadata, URLs, README, and NEWS for the GitHub workflow;
* adds fixed-seed numerical regression tests for core estimation, regularization,
  derivative estimation, and bandwidth selection;
* improves selected internal matrix and covariance computations while preserving
  the numerical regression baseline;
* regenerates and cleans package documentation.

This is a resubmission. In the first submission, CRAN incoming checks reported
non-canonical or permanently moved URLs in README.md. These URLs have been
updated to use the canonical CRAN package page and current NSF award URLs.

The previous CRAN maintainer was Rajita Chandak. The current maintainer is
Matias D. Cattaneo.

## R CMD check results

Checked locally on Windows 11 x64 using R version 4.6.0 (2026-04-24 ucrt).

0 errors | 0 warnings | 2 notes

The first NOTE is the CRAN incoming feasibility NOTE for the maintainer change:

* New maintainer: Matias D. Cattaneo <matias.d.cattaneo@gmail.com>
* Old maintainer: Rajita Chandak <rajita.chandak@epfl.ch>

The second NOTE is local to this check environment:

* Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being
  installed.

Pandoc is not available on this local Windows machine; the package README and
NEWS files are plain Markdown files and are expected to be checked normally on
CRAN's check systems.

## Reverse dependencies

There are currently no downstream dependencies for this package.
