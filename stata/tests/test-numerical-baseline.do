********************************************************************************
* LPCDE STATA PACKAGE -- Numerical baseline smoke test
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

quietly lpcde y x, grid(ygrid) xeval(0) bw(1)
matrix R = e(result)

scalar tol = 1e-7
matrix E = (-2, 1, -0.035998710806018033, -0.0055464210156906078, 0.28614641373295102, 0.31298560947401322 \ ///
            -1, 1,  0.34450258455391702,   0.41330101026916793,   0.035701793707458795, 0.083307920859939932 \ ///
             0, 1,  0.32968206975370806,   0.39711767972091327,   0.016539301915963513, 0.064233809916422191 \ ///
             1, 1,  0.22507194337292116,   0.16731018717777676,   0.025224524918620032, 0.070304640197867016 \ ///
             2, 1,  0.073006038706733986,  0.06702511721084832,   0.050384629913758666, 0.15654752908295377)

forvalues i = 1/5 {
    assert abs(el(R, `i', 1) - el(E, `i', 1)) < tol
    assert abs(el(R, `i', 2) - el(E, `i', 2)) < tol
    assert abs(el(R, `i', 4) - el(E, `i', 3)) < tol
    assert abs(el(R, `i', 5) - el(E, `i', 4)) < tol
    assert abs(el(R, `i', 6) - el(E, `i', 5)) < tol
    assert abs(el(R, `i', 7) - el(E, `i', 6)) < tol
}

quietly lpbwcde y x, grid(ygrid) xeval(0)
matrix B = e(result)
matrix EB = (-2, 1.2319871443245998, 25 \ ///
             -1, 1.2319871443245998, 95 \ ///
              0, 1.2319871443245998, 137 \ ///
              1, 1.2319871443245998, 92 \ ///
              2, 1.2319871443245998, 29)
forvalues i = 1/5 {
    assert abs(el(B, `i', 1) - el(EB, `i', 1)) < tol
    assert abs(el(B, `i', 2) - el(EB, `i', 2)) < tol
    assert abs(el(B, `i', 3) - el(EB, `i', 3)) < tol
}

display as text "Standalone Stata lpcde numerical smoke fixture OK"
