********************************************************************************
* LPCDE STATA PACKAGE -- lpbwcde
* Authors: Matias D. Cattaneo, Rajita Chandak, Michael Jansson, Xinwei Ma
********************************************************************************
*!version 1.0.0 2026-05-21

capture program drop lpbwcde

program define lpbwcde, eclass
	version 16.0
	syntax varlist(min=2 numeric) [if] [in] [, ///
		GRId(varname numeric) ///
		XEval(numlist) ///
		P(integer -1) ///
		Q(integer -1) ///
		MU(integer -1) ///
		NU(integer -1) ///
		KERnel(string) ///
		BWSelect(string) ///
		NG(integer -1) ///
		GRIDSpacing(string) ///
		noREGularize ///
		GENvars(name) ///
		RGRId(varname numeric) ///
		RINDex(varname numeric) ///
		SEParator(integer 5) ///
		]

	marksample touse
	gettoken yvar xvars : varlist
	local dx : word count `xvars'

	qui count if `touse'
	local n = r(N)
	if (`n' == 0) {
		di as err "no observations"
		exit 2000
	}

	if (`mu' < 0) local mu_default = 1
	else local mu_default = 0
	if (`nu' < 0) local nu_default = 1
	else local nu_default = 0
	if (`p' < 0) {
		if (`mu_default') local p = 2
		else local p = `mu' + 1
	}
	if (`q' < 0) {
		if (`nu_default') local q = 1
		else local q = `nu' + 1
	}
	if (`mu' < 0) local mu = min(1, `p')
	if (`nu' < 0) local nu = min(0, `q')

	foreach opt in p q mu nu {
		if (``opt'' < 0 | ``opt'' > 20) {
			di as err "`opt'(): has to be an integer between 0 and 20"
			exit 198
		}
	}
	if (`mu' > `p' | `nu' > `q') {
		di as err "mu() and nu() cannot exceed p() and q()"
		exit 198
	}
	if (`dx' > 1 & `q' >= 5) {
		di as err "multivariate conditioning polynomials of order 5 or larger are not supported"
		exit 198
	}

	if ("`kernel'" == "") local kernel = "epanechnikov"
	else {
		local kernel = lower("`kernel'")
		if ("`kernel'" != "epanechnikov" & "`kernel'" != "triangular" & "`kernel'" != "uniform") {
			di as err "kernel(): incorrectly specified: options(epanechnikov, triangular, uniform)"
			exit 198
		}
	}

	if ("`bwselect'" == "") local bwselect = "imse-rot"
	else {
		local bwselect = lower("`bwselect'")
		if ("`bwselect'" != "mse-rot" & "`bwselect'" != "imse-rot") {
			di as err "bwselect(): incorrectly specified: options(mse-rot, imse-rot)"
			exit 198
		}
	}

	local has_xeval = 0
	if ("`xeval'" != "") {
		local nxeval : word count `xeval'
		if (`nxeval' != `dx') {
			di as err "xeval(): must contain one number for each conditioning variable"
			exit 198
		}
		matrix Xeval = (`xeval')
		local has_xeval = 1
	}
	if (`separator' <= 1) local separator = 1
	local regularizeflag = cond("`regularize'" == "", 1, 0)

	if ("`genvars'" != "") {
		foreach suffix in grid bw nh {
			confirm new variable `genvars'_`suffix'
		}
	}

	capture findfile lpcde_fun.do
	if (_rc) {
		di as err "lpcde_fun.do not found"
		exit 601
	}
	quietly do "`r(fn)'"

	tempvar temp_touse
	qui gen byte `temp_touse' = `touse'

	mata: lpcde_stata_bw_run("`varlist'", "`temp_touse'", "`grid'", "`rgrid'", "`rindex'", `has_xeval', `ng', "`gridspacing'", `p', `q', `mu', `nu', "`kernel'", "`bwselect'", `regularizeflag')

	matrix colnames Result = y_grid bw nh index
	local ng_result = rowsof(Result)

	if ("`genvars'" != "") {
		foreach suffix in grid bw nh {
			qui gen double `genvars'_`suffix' = .
		}
		forvalues i = 1/`ng_result' {
			qui replace `genvars'_grid = Result[`i', 1] in `i'
			qui replace `genvars'_bw = Result[`i', 2] in `i'
			qui replace `genvars'_nh = Result[`i', 3] in `i'
		}
		label variable `genvars'_grid "lpbwcde: y grid point"
		label variable `genvars'_bw "lpbwcde: bandwidth"
		label variable `genvars'_nh "lpbwcde: effective sample size"
	}

	disp ""
	disp "Bandwidth Selection for Local Polynomial Conditional Density Estimation."
	disp ""
	disp in smcl in gr "{lalign 1: Sample size                              (n=)    }" _col(19) in ye %15.0f `n'
	disp in smcl in gr "{lalign 1: Conditioning variables                   (d=)    }" _col(19) in ye %15.0f `dx'
	disp in smcl in gr "{lalign 1: Polynomial order in y                    (p=)    }" _col(19) in ye %15.0f `p'
	disp in smcl in gr "{lalign 1: Polynomial order in x                    (q=)    }" _col(19) in ye %15.0f `q'
	disp in smcl in gr "{lalign 1: Derivative order in y                   (mu=)    }" _col(19) in ye %15.0f `mu'
	disp in smcl in gr "{lalign 1: Derivative order in x                   (nu=)    }" _col(19) in ye %15.0f `nu'
	disp in smcl in gr "{lalign 1: Kernel function                                  }" _col(19) in ye "{ralign 15: `kernel'}"
	disp in smcl in gr "{lalign 1: Bandwidth selection method                       }" _col(19) in ye "{ralign 15: `bwselect'}"
	disp ""

	disp in smcl in gr "{hline 42}"
	disp in smcl in gr "{ralign 4:Index}" _col(7) "{ralign 9:y_grid}" _col(17) "{ralign 9:B.W.}" _col(27) "{ralign 8:Eff.n}"
	disp in smcl in gr "{hline 42}"
	forvalues i = 1/`ng_result' {
		disp in smcl in gr %4.0f Result[`i', 4] _col(7) in ye %9.4f Result[`i', 1] _col(17) %9.4f Result[`i', 2] _col(27) %8.0f Result[`i', 3]
		if (`i' != `ng_result' & `separator' > 1 & mod(`i', `separator') == 0) {
			disp in smcl in gr "{hline 42}"
		}
	}
	disp in smcl in gr "{hline 42}"

	ereturn clear
	ereturn scalar N = `n'
	ereturn scalar d = `dx'
	ereturn scalar p = `p'
	ereturn scalar q = `q'
	ereturn scalar mu = `mu'
	ereturn scalar nu = `nu'
	ereturn local yvar = "`yvar'"
	ereturn local xvars = "`xvars'"
	ereturn local kernel = "`kernel'"
	ereturn local bwselect = "`bwselect'"
	ereturn matrix result = Result
end
