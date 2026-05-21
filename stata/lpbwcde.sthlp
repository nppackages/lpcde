{smcl}
{* *!version 1.0.0 2026-05-21}{...}

{title:Title}

{p 4 8}{cmd:lpbwcde} {hline 2} Bandwidth selection for local polynomial conditional density estimation.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:lpbwcde} {it:yvar xvars} {ifin}
[{cmd:,}
{cmd:grid(}{it:var}{cmd:)}
{cmd:xeval(}{it:numlist}{cmd:)}
{cmd:p(}{it:#}{cmd:)}
{cmd:q(}{it:#}{cmd:)}
{cmd:mu(}{it:#}{cmd:)}
{cmd:nu(}{it:#}{cmd:)}
{cmd:kernel(}{it:KernelFn}{cmd:)}
{cmd:bwselect(}{it:BwMethod}{cmd:)}
{cmd:ng(}{it:#}{cmd:)}
{cmd:gridspacing(}{it:Spacing}{cmd:)}
{cmd:noregularize}
{cmd:genvars(}{it:prefix}{cmd:)}
{cmd:rgrid(}{it:var}{cmd:)}
{cmd:rindex(}{it:var}{cmd:)}
{cmd:separator(}{it:#}{cmd:)}
]{p_end}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:lpbwcde} implements rule-of-thumb bandwidth selection for local polynomial conditional density estimation.
The first variable is the outcome and all remaining variables are conditioning covariates. Numerical calculations are implemented in standalone Mata code.{p_end}

{p 8 8}Companion command: {help lpcde:lpcde} for estimation and inference.{p_end}

{marker options}{...}
{title:Options}

{p 4 8}{opt grid(var)} specifies the y-grid for bandwidth selection. When omitted, the R package default grid is used.{p_end}

{p 4 8}{opt xeval(numlist)} specifies the conditioning value. The list must contain one number for each conditioning covariate. When omitted, componentwise medians are used.{p_end}

{p 4 8}{opt p(#)} and {opt q(#)} specify local polynomial orders in the y- and x-directions. Defaults are {cmd:p(2)} and {cmd:q(1)}.{p_end}

{p 4 8}{opt mu(#)} and {opt nu(#)} specify derivative orders with respect to y and x. Defaults are {cmd:mu(1)} and {cmd:nu(0)} subject to the polynomial orders.{p_end}

{p 4 8}{opt kernel(KernelFn)} specifies {cmd:epanechnikov}, {cmd:triangular}, or {cmd:uniform}. Default is {cmd:epanechnikov}.{p_end}

{p 4 8}{opt bwselect(BwMethod)} specifies {cmd:imse-rot} or {cmd:mse-rot}. Default is {cmd:imse-rot}.{p_end}

{p 4 8}{opt ng(#)} and {opt gridspacing(Spacing)} control automatic grid generation when {cmd:grid()} is omitted. {cmd:gridspacing(quantile)} requests quantile spacing.{p_end}

{p 4 8}{opt noregularize} suppresses the local sample size regularization used by the default bandwidth selector.{p_end}

{p 4 8}{opt genvars(prefix)} stores double-precision results in variables named {it:prefix}{cmd:_grid}, {it:prefix}{cmd:_bw}, and {it:prefix}{cmd:_nh}.{p_end}

{p 4 8}{opt rgrid(var)} displays rows nearest to the supplied y-grid. {opt rindex(var)} displays rows nearest to supplied row indices.{p_end}

{p 4 8}{opt separator(#)} draws a separator line after every {it:#} rows in the display.{p_end}

{marker saved_results}{...}
{title:Saved Results}

{p 4 8}{cmd:lpbwcde} saves the following in {cmd:e()}:{p_end}

{synoptset 20 tabbed}{...}
{synopt:{cmd:e(result)}}matrix with columns {cmd:y_grid bw nh index}{p_end}
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
{p 4 8}{cmd:. lpbwcde y x, grid(ygrid) xeval(0) genvars(bw)}{p_end}

{title:Authors}

{p 4 8}Matias D. Cattaneo, Princeton University.
{browse "mailto:matias.d.cattaneo@gmail.com":matias.d.cattaneo@gmail.com}.{p_end}

{p 4 8}Rajita Chandak, University of Wisconsin-Madison.
{browse "mailto:rajita.chandak@gmail.com":rajita.chandak@gmail.com}.{p_end}

{p 4 8}Michael Jansson, University of California, Berkeley.
{browse "mailto:michael.jansson.berkeley@gmail.com":michael.jansson.berkeley@gmail.com}.{p_end}

{p 4 8}Xinwei Ma, University of California, San Diego.
{browse "mailto:xinweima.pku@gmail.com":xinweima.pku@gmail.com}.{p_end}
