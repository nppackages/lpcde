********************************************************************************
* LPCDE Stata replication script for the software article examples
********************************************************************************
version 16.0
clear all
set more off

adopath ++ "stata"

import delimited "Python/lpcde/tests/fixtures/r-baseline-1d.csv", clear
recast double x y
generate double ygrid = .
forvalues i = 1/9 {
    replace ygrid = -2 + (`i' - 1) * 0.5 in `i'
}

lpcde y x, grid(ygrid) xeval(0) bw(1)
matrix list e(result)

lpcde y x, grid(ygrid) xeval(0) bw(1) nonneg normalize
matrix list e(result)

lpcde y x, grid(ygrid) xeval(0) bw(1) p(3) mu(2)
matrix list e(result)

lpbwcde y x, grid(ygrid) xeval(0)
matrix list e(result)
