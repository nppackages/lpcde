********************************************************************************
* LPCDE STATA PACKAGE -- lpcde
* Authors: Matias D. Cattaneo, Rajita Chandak, Michael Jansson, Xinwei Ma
********************************************************************************
*!version 1.0.0 2026-05-21

capture program drop lpcde

program define lpcde, eclass
	version 16.0
	syntax varlist(min=2 numeric) [if] [in] [, ///
		GRId(varname numeric) ///
		XEval(numlist) ///
		BW(string) ///
		P(integer -1) ///
		Q(integer -1) ///
		PRBC(integer -1) ///
		QRBC(integer -1) ///
		MU(integer -1) ///
		NU(integer -1) ///
		noRBC ///
		KERnel(string) ///
		BWSelect(string) ///
		COVFlag(string) ///
		NG(integer -1) ///
		GRIDSpacing(string) ///
		NONNEG ///
		NORMalize ///
		GENvars(name) ///
		RGRId(varname numeric) ///
		RINDex(varname numeric) ///
		Level(real 95) ///
		SEParator(integer 5) ///
		PLot ///
		GRAph_options(string) ///
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
	if (`prbc' < 0) local prbc = `p' + 1
	if (`qrbc' < 0) local qrbc = `q' + 1

	foreach opt in p q prbc qrbc mu nu {
		if (``opt'' < 0 | ``opt'' > 20) {
			di as err "`opt'(): has to be an integer between 0 and 20"
			exit 198
		}
	}
	if (`prbc' < `p' | `qrbc' < `q') {
		di as err "prbc() and qrbc() must be at least as large as p() and q()"
		exit 198
	}
	if (`mu' > `p' | `nu' > `q') {
		di as err "mu() and nu() cannot exceed p() and q()"
		exit 198
	}
	if (`dx' > 1 & max(`q', `qrbc') >= 5) {
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

	if ("`covflag'" == "") local covflag = "full"
	else {
		local covflag = lower("`covflag'")
		if ("`covflag'" != "full" & "`covflag'" != "diag" & "`covflag'" != "off") {
			di as err "covflag(): incorrectly specified: options(full, diag, off)"
			exit 198
		}
	}

	if (`level' <= 0 | `level' >= 100) {
		di as err "level(): incorrectly specified: should be between 0 and 100"
		exit 198
	}
	if (`separator' <= 1) local separator = 1

	local rbcflag = cond("`rbc'" == "", 1, 0)
	local nonnegflag = cond("`nonneg'" == "", 0, 1)
	local normalizeflag = cond("`normalize'" == "", 0, 1)

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

	local has_bw = 0
	local bw_is_var = 0
	local bw_scalar = .
	if ("`bw'" != "") {
		local has_bw = 1
		capture confirm number `bw'
		if (_rc == 0) {
			if (`bw' <= 0) {
				di as err "bw(): must be positive"
				exit 198
			}
			local bw_scalar = `bw'
		}
		else {
			capture confirm numeric variable `bw'
			if (_rc) {
				di as err "bw(): must be a positive scalar or numeric variable"
				exit 198
			}
			local bw_is_var = 1
		}
	}

	if ("`genvars'" != "") {
		foreach suffix in grid bw nh est est_rbc se se_rbc ci_l ci_r {
			confirm new variable `genvars'_`suffix'
		}
	}

	capture mata: mata which lpcde_stata_run()
	if (_rc) {
		capture quietly mata: mata mlib index
		capture mata: mata which lpcde_stata_run()
	}
	if (_rc) {
		di as err "lpcde.mlib not found or not indexed"
		di as err `"rebuild with Stata 16: do "stata/build-lpcde-mlib.do""'
		exit 601
	}

	tempvar temp_touse
	qui gen byte `temp_touse' = `touse'

	mata: lpcde_stata_run("`varlist'", "`temp_touse'", "`grid'", "`rgrid'", "`rindex'", "`bw'", `has_bw', `bw_is_var', `bw_scalar', `has_xeval', `ng', "`gridspacing'", `p', `q', `prbc', `qrbc', `mu', `nu', `rbcflag', "`kernel'", "`bwselect'", "`covflag'", `nonnegflag', `normalizeflag', `level')

	matrix colnames Result = y_grid bw nh est est_RBC se se_RBC CI_l CI_r index
	local ng_result = rowsof(Result)

	if ("`genvars'" != "") {
		foreach suffix in grid bw nh est est_rbc se se_rbc ci_l ci_r {
			qui gen double `genvars'_`suffix' = .
		}
		forvalues i = 1/`ng_result' {
			qui replace `genvars'_grid = Result[`i', 1] in `i'
			qui replace `genvars'_bw = Result[`i', 2] in `i'
			qui replace `genvars'_nh = Result[`i', 3] in `i'
			qui replace `genvars'_est = Result[`i', 4] in `i'
			qui replace `genvars'_est_rbc = Result[`i', 5] in `i'
			qui replace `genvars'_se = Result[`i', 6] in `i'
			qui replace `genvars'_se_rbc = Result[`i', 7] in `i'
			qui replace `genvars'_ci_l = Result[`i', 8] in `i'
			qui replace `genvars'_ci_r = Result[`i', 9] in `i'
		}
		label variable `genvars'_grid "lpcde: y grid point"
		label variable `genvars'_bw "lpcde: bandwidth"
		label variable `genvars'_nh "lpcde: effective sample size"
		label variable `genvars'_est "lpcde: point estimate"
		label variable `genvars'_est_rbc "lpcde: robust bias-corrected estimate"
		label variable `genvars'_se "lpcde: standard error"
		label variable `genvars'_se_rbc "lpcde: robust bias-corrected standard error"
		label variable `genvars'_ci_l "lpcde: left level-`level' confidence interval"
		label variable `genvars'_ci_r "lpcde: right level-`level' confidence interval"
	}

	disp ""
	disp "Local Polynomial Conditional Density Estimation and Inference."
	disp ""
	disp in smcl in gr "{lalign 1: Sample size                              (n=)    }" _col(19) in ye %15.0f `n'
	disp in smcl in gr "{lalign 1: Conditioning variables                   (d=)    }" _col(19) in ye %15.0f `dx'
	disp in smcl in gr "{lalign 1: Polynomial order in y                    (p=)    }" _col(19) in ye %15.0f `p'
	disp in smcl in gr "{lalign 1: Polynomial order in x                    (q=)    }" _col(19) in ye %15.0f `q'
	disp in smcl in gr "{lalign 1: Derivative order in y                   (mu=)    }" _col(19) in ye %15.0f `mu'
	disp in smcl in gr "{lalign 1: Derivative order in x                   (nu=)    }" _col(19) in ye %15.0f `nu'
	disp in smcl in gr "{lalign 1: Kernel function                                  }" _col(19) in ye "{ralign 15: `kernel'}"
	disp in smcl in gr "{lalign 1: Bandwidth selection method                       }" _col(19) in ye "{ralign 15: `bwselect'}"
	disp in smcl in gr "{lalign 1: Covariance computation                           }" _col(19) in ye "{ralign 15: `covflag'}"
	disp ""

	disp in smcl in gr "{hline 82}"
	disp in smcl in gr "{ralign 4:Index}" _col(7) "{ralign 9:y_grid}" _col(17) "{ralign 9:B.W.}" _col(27) "{ralign 7:Eff.n}" _col(36) "{ralign 10:Est.}" _col(48) "{ralign 10:Std.Err.}" _col(60) "{ralign 20:`level'% C.I.}"
	disp in smcl in gr "{hline 82}"
	forvalues i = 1/`ng_result' {
		disp in smcl in gr %4.0f Result[`i', 10] _col(7) in ye %9.4f Result[`i', 1] _col(17) %9.4f Result[`i', 2] _col(27) %7.0f Result[`i', 3] _col(36) %10.4f Result[`i', 4] _col(48) %10.4f Result[`i', 6] _col(60) %10.4f Result[`i', 8] _col(71) %10.4f Result[`i', 9]
		if (`i' != `ng_result' & `separator' > 1 & mod(`i', `separator') == 0) {
			disp in smcl in gr "{hline 82}"
		}
	}
	disp in smcl in gr "{hline 82}"

	if ("`plot'" != "") {
		tempvar plot_grid plot_est plot_cil plot_cir
		qui gen double `plot_grid' = .
		qui gen double `plot_est' = .
		qui gen double `plot_cil' = .
		qui gen double `plot_cir' = .
		forvalues i = 1/`ng_result' {
			qui replace `plot_grid' = Result[`i', 1] in `i'
			qui replace `plot_est' = Result[`i', 4] in `i'
			qui replace `plot_cil' = Result[`i', 8] in `i'
			qui replace `plot_cir' = Result[`i', 9] in `i'
		}
		if (`"`graph_options'"' == "") {
			local graph_options = `"legend(off) title("lpcde (p=`p', q=`q', mu=`mu', nu=`nu')") xtitle("y grid") ytitle("")"'
		}
		twoway (rarea `plot_cil' `plot_cir' `plot_grid', sort color(red%30)) ///
			(line `plot_est' `plot_grid', sort lcolor(red) lwidth(medthin)), ///
			`graph_options'
	}

	ereturn clear
	ereturn scalar N = `n'
	ereturn scalar d = `dx'
	ereturn scalar p = `p'
	ereturn scalar q = `q'
	ereturn scalar p_RBC = `prbc'
	ereturn scalar q_RBC = `qrbc'
	ereturn scalar mu = `mu'
	ereturn scalar nu = `nu'
	ereturn scalar level = `level'
	ereturn local yvar = "`yvar'"
	ereturn local xvars = "`xvars'"
	ereturn local kernel = "`kernel'"
	ereturn local bwselect = "`bwselect'"
	ereturn local covflag = "`covflag'"
	ereturn matrix result = Result
	ereturn matrix CovMat = CovMat
	ereturn matrix CovMat_RBC = CovMat_RBC
end
