# Changelog

Notable project changes are listed from newest to oldest.

## 2026-05-21 - Version 1.0.0

- Added a GitHub Actions workflow for PyPI Trusted Publishing of the Python package on published releases.
- Normalized author affiliations to university names only across package metadata and Stata help files.
- Added Stata help-file author sections matching the LPDENSITY help-file style and removed the dependency section from the Stata help surface.
- Added Stata help PDFs for `lpcde` and `lpbwcde`, generated from the `.sthlp` files using the LPDENSITY-style local build script.
- Updated the Stata illustration do-file so it can run directly from a local repository checkout.
- Simplified the Stata API by removing the `precision()` option; generated variables are now always double precision.
- Renamed `replication_scripts` to `replication` and added Python and Stata analogues for the illustration and software article replication entry points.
- Bumped the R package and Python package versions to `1.0.0`.
- Marked this checkpoint as the first cross-language GitHub modernization release, with R numerical baselines and Python/R alignment tests in place.
- Added the initial standalone Stata/Mata package under `stata`, including `lpcde`/`lpbwcde` commands, help files, package metadata, an illustration do-file, and a numerical smoke baseline for cross-language parity.

## 2026-05-21 - Python Package Scaffold And Alignment

- Added the initial Python package under `Python/lpcde`, including package metadata, README, and public `lpcde()`/`lpbwcde()` APIs.
- Ported the R estimator, covariance, bandwidth-selection, and result-summary paths to Python using NumPy, SciPy, and pandas.
- Restructured the Python package to match the `lpdensity` layout, with command-level modules, function-helper modules, validation helpers, requirements file, and an illustration script under `Python/`.
- Added R-generated Python test fixtures and a pytest numerical alignment suite covering fixed bandwidth estimation, regularized output, derivatives, kernel variants, covariance flags, two-dimensional conditioning, and bandwidth selectors.
- Extended CI and repository layout checks so Python package tests run alongside the R package checks.

## 2026-05-21 - R Package Cleanup And Speed Baseline

- Refreshed package help generation and fixed the `coef.lpbwcde()` examples so they are emitted as examples rather than malformed `seealso` content.
- Cleaned R package metadata by removing non-portable compiler overrides from `src/Makevars` and `src/Makevars.win`, and by dropping unused DESCRIPTION fields.
- Updated the package-level README and NEWS entries to match the GitHub-oriented repository layout and current package metadata.
- Optimized selected one-dimensional matrix and covariance paths while preserving the fixed-seed numerical regression baseline.
- Revalidated the R package with numerical baseline checks and local package build/check commands.

## 2026-05-20 - GitHub Workflow Setup

- Reorganized the top-level README to match the current `binsreg` layout, with package command summaries, R installation notes, replication links, references, and funding.
- Standardized repository maintenance files for GitHub use, including CI, issue templates, security policy, pull request template, and Dependabot configuration.
- Updated R package ownership metadata so Matias D. Cattaneo is the maintainer and contributor email addresses match the current project records.
- Added a fixed-seed R numerical regression fixture and testthat check to protect core `lpcde`, regularized `lpcde`, derivative, and `lpbwcde` outputs from unintentional drift.
- Standardized the project license file with the short GPL-3 license notice used by `binsreg`.
