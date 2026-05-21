{smcl}
{* *!version 1.0.0 2026-05-21}{...}

{title:Title}

{p 4 8}{cmd:lpcde} {hline 2} Local polynomial conditional density estimation and inference.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:lpcde} {it:yvar xvars} {ifin}
[{cmd:,}
{cmd:grid(}{it:var}{cmd:)}
{cmd:xeval(}{it:numlist}{cmd:)}
{cmd:bw(}{it:var} or {it:#}{cmd:)}
{cmd:p(}{it:#}{cmd:)}
{cmd:q(}{it:#}{cmd:)}
{cmd:prbc(}{it:#}{cmd:)}
{cmd:qrbc(}{it:#}{cmd:)}
{cmd:mu(}{it:#}{cmd:)}
{cmd:nu(}{it:#}{cmd:)}
{cmd:norbc}
{cmd:kernel(}{it:KernelFn}{cmd:)}
{cmd:bwselect(}{it:BwMethod}{cmd:)}
{cmd:covflag(}{it:CovMethod}{cmd:)}
{cmd:ng(}{it:#}{cmd:)}
{cmd:gridspacing(}{it:Spacing}{cmd:)}
{cmd:nonneg}
{cmd:normalize}
{cmd:genvars(}{it:prefix}{cmd:)}
{cmd:rgrid(}{it:var}{cmd:)}
{cmd:rindex(}{it:var}{cmd:)}
{cmd:level(}{it:#}{cmd:)}
{cmd:separator(}{it:#}{cmd:)}
{cmd:plot}
{cmd:graph_options(}{it:GraphOpts}{cmd:)}
]{p_end}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:lpcde} implements local polynomial conditional CDF, density, and derivative estimation.
The first variable in {it:yvar xvars} is the outcome; remaining variables are the conditioning covariates.
Numerical calculations are implemented in standalone Mata code and are intended to match the companion R and Python packages.{p_end}

{p 8 8}Companion command: {help lpbwcde:lpbwcde} for bandwidth selection.{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Estimation}

{p 4 8}{opt grid(var)} specifies the y-grid on which the conditional density is estimated. When omitted, the R package default grid is used.{p_end}

{p 4 8}{opt xeval(numlist)} specifies the conditioning value. The list must contain one number for each conditioning covariate. When omitted, componentwise medians are used.{p_end}

{p 4 8}{opt bw(var|#)} specifies a bandwidth variable or a positive scalar bandwidth. When omitted, {help lpbwcde:lpbwcde} is called internally with {cmd:bwselect()}.{p_end}

{p 4 8}{opt p(#)} and {opt q(#)} specify local polynomial orders in the y- and x-directions. Defaults are {cmd:p(2)} and {cmd:q(1)}.{p_end}

{p 4 8}{opt prbc(#)} and {opt qrbc(#)} specify robust bias-correction polynomial orders. Defaults are {cmd:p()+1} and {cmd:q()+1}.{p_end}

{p 4 8}{opt mu(#)} and {opt nu(#)} specify derivative orders with respect to y and x. Defaults are {cmd:mu(1)} and {cmd:nu(0)} subject to the polynomial orders.{p_end}

{p 4 8}{opt norbc} suppresses robust bias correction and centers confidence intervals at the conventional estimate.{p_end}

{p 4 8}{opt kernel(KernelFn)} specifies {cmd:epanechnikov}, {cmd:triangular}, or {cmd:uniform}. Default is {cmd:epanechnikov}.{p_end}

{p 4 8}{opt covflag(CovMethod)} specifies covariance computation: {cmd:full}, {cmd:diag}, or {cmd:off}. Default is {cmd:full}.{p_end}

{p 4 8}{opt nonneg} truncates conventional point estimates below at zero. {opt normalize} normalizes conventional estimates to integrate to one over the y-grid.{p_end}

{dlgtab:Bandwidth Selection}

{p 4 8}{opt bwselect(BwMethod)} specifies {cmd:imse-rot} or {cmd:mse-rot}. Default is {cmd:imse-rot}; ignored when {cmd:bw()} is supplied.{p_end}

{p 4 8}{opt ng(#)} and {opt gridspacing(Spacing)} control automatic grid generation when {cmd:grid()} is omitted. {cmd:gridspacing(quantile)} requests quantile spacing.{p_end}

{dlgtab:Storing and Displaying Results}

{p 4 8}{opt genvars(prefix)} stores double-precision results in variables named {it:prefix}{cmd:_grid}, {it:prefix}{cmd:_bw}, {it:prefix}{cmd:_nh}, {it:prefix}{cmd:_est}, {it:prefix}{cmd:_est_rbc}, {it:prefix}{cmd:_se}, {it:prefix}{cmd:_se_rbc}, {it:prefix}{cmd:_ci_l}, and {it:prefix}{cmd:_ci_r}.{p_end}

{p 4 8}{opt rgrid(var)} displays rows nearest to the supplied y-grid. {opt rindex(var)} displays rows nearest to supplied row indices.{p_end}

{p 4 8}{opt level(#)} controls pointwise confidence interval level. Default is {cmd:level(95)}.{p_end}

{p 4 8}{opt separator(#)} draws a separator line after every {it:#} rows in the display.{p_end}

{dlgtab:Plotting}

{p 4 8}{opt plot} plots the conventional point estimate with pointwise confidence intervals. {opt graph_options()} passes options to {cmd:twoway}.{p_end}

{marker saved_results}{...}
{title:Saved Results}

{p 4 8}{cmd:lpcde} saves the following in {cmd:e()}:{p_end}

{synoptset 20 tabbed}{...}
{synopt:{cmd:e(result)}}matrix with columns {cmd:y_grid bw nh est est_RBC se se_RBC CI_l CI_r index}{p_end}
{synopt:{cmd:e(CovMat)}}conventional covariance matrix{p_end}
{synopt:{cmd:e(CovMat_RBC)}}robust bias-corrected covariance matrix{p_end}
{synopt:{cmd:e(N)}}sample size{p_end}
{synopt:{cmd:e(d)}}number of conditioning covariates{p_end}
{synopt:{cmd:e(p)}, {cmd:e(q)}, {cmd:e(mu)}, {cmd:e(nu)}}polynomial and derivative orders{p_end}

{marker examples}{...}
{title:Examples}

{p 4 8}{cmd:. set obs 500}{p_end}
{p 4 8}{cmd:. set seed 42}{p_end}
{p 4 8}{cmd:. gen double x = rnormal()}{p_end}
{p 4 8}{cmd:. gen double y = x + rnormal()}{p_end}
{p 4 8}{cmd:. gen double ygrid = -2 + (_n-1) if _n <= 5}{p_end}
{p 4 8}{cmd:. lpcde y x, grid(ygrid) xeval(0) bw(1) genvars(cde)}{p_end}

{title:Authors}

{p 4 8}Matias D. Cattaneo, Princeton University.
{browse "mailto:matias.d.cattaneo@gmail.com":matias.d.cattaneo@gmail.com}.{p_end}

{p 4 8}Rajita Chandak, University of Wisconsin-Madison.
{browse "mailto:rajita.chandak@gmail.com":rajita.chandak@gmail.com}.{p_end}

{p 4 8}Michael Jansson, University of California, Berkeley.
{browse "mailto:michael.jansson.berkeley@gmail.com":michael.jansson.berkeley@gmail.com}.{p_end}

{p 4 8}Xinwei Ma, University of California, San Diego.
{browse "mailto:xinweima.pku@gmail.com":xinweima.pku@gmail.com}.{p_end}
