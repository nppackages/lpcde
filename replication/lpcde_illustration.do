********************************************************************************
* LPCDE Stata replication illustration
********************************************************************************
version 16.0
clear all
set more off

adopath ++ "stata"

import delimited "Python/lpcde/tests/fixtures/r-baseline-1d.csv", clear
recast double x y
generate double ygrid = .
replace ygrid = -2 in 1
replace ygrid = -1 in 2
replace ygrid = 0 in 3
replace ygrid = 1 in 4
replace ygrid = 2 in 5

lpcde y x, grid(ygrid) xeval(0) bw(1) genvars(cde)
matrix list e(result)

lpbwcde y x, grid(ygrid) xeval(0) genvars(bw)
matrix list e(result)
