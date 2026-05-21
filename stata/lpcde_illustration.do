********************************************************************************
* LPCDE STATA PACKAGE -- Illustration
********************************************************************************
version 16.0
clear all
set more off
set seed 42

local pwd = c(pwd)
capture confirm file "`pwd'/stata/lpcde.ado"
if !_rc {
	adopath ++ "`pwd'/stata"
}
else {
	capture confirm file "`pwd'/lpcde.ado"
	if !_rc {
		adopath ++ "`pwd'"
	}
}

set obs 500
generate double x = rnormal()
generate double y = x + rnormal()
generate double ygrid = -2 + (_n - 1) if _n <= 5

lpbwcde y x, grid(ygrid) xeval(0) genvars(bw)
lpcde y x, grid(ygrid) xeval(0) bw(1) genvars(cde)

matrix list e(result)
