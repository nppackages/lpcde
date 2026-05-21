********************************************************************************
* LPCDE STATA PACKAGE -- standalone Mata functions
* Authors: Matias D. Cattaneo, Rajita Chandak, Michael Jansson, Xinwei Ma
********************************************************************************
*!version 1.0.0 2026-05-21

capture mata mata drop lpcde_*()

mata

real scalar lpcde_factorial(real scalar n)
{
	real scalar out, i
	out = 1
	if (n <= 1) return(out)
	for (i=2; i<=n; i++) out = out * i
	return(out)
}

real scalar lpcde_sd(real colvector x)
{
	real scalar n
	n = rows(x)
	if (n <= 1) return(0)
	return(sqrt(sum((x :- mean(x)):^2) / (n - 1)))
}

real rowvector lpcde_colsd(real matrix X)
{
	real scalar j
	real rowvector out
	out = J(1, cols(X), .)
	for (j=1; j<=cols(X); j++) out[j] = lpcde_sd(X[., j])
	return(out)
}

real rowvector lpcde_colmean(real matrix X)
{
	real scalar j
	real rowvector out
	out = J(1, cols(X), .)
	for (j=1; j<=cols(X); j++) out[j] = mean(X[., j])
	return(out)
}

real scalar lpcde_rowmean(real rowvector x)
{
	return(sum(x) / cols(x))
}

real matrix lpcde_sweep_row(real matrix X, real rowvector v, string scalar op)
{
	real matrix out
	real scalar j
	out = X
	for (j=1; j<=cols(X); j++) {
		if (op == "-") out[., j] = X[., j] :- v[j]
		else if (op == "/") out[., j] = X[., j] :/ v[j]
		else if (op == "*") out[., j] = X[., j] :* v[j]
	}
	return(out)
}

real colvector lpcde_nonmissing(real colvector x)
{
	return(select(x, x :< .))
}

real colvector lpcde_which(real colvector mask)
{
	return(select((1::rows(mask)), mask))
}

real colvector lpcde_intersect(real colvector a, real colvector b)
{
	real scalar i, j, n
	real colvector out
	out = J(rows(a), 1, .)
	n = 0
	for (i=1; i<=rows(a); i++) {
		for (j=1; j<=rows(b); j++) {
			if (a[i] == b[j]) {
				n++
				out[n] = a[i]
				break
			}
		}
	}
	if (n == 0) return(J(0, 1, .))
	return(out[1..n])
}

real scalar lpcde_quantile(real colvector x, real scalar p)
{
	real scalar n, pos, lo, hi, w
	real colvector s
	s = sort(x, 1)
	n = rows(s)
	if (n == 1) return(s[1])
	pos = 1 + (n - 1) * p
	lo = floor(pos)
	hi = ceil(pos)
	w = pos - lo
	if (lo == hi) return(s[lo])
	return((1 - w) * s[lo] + w * s[hi])
}

real colvector lpcde_linspace(real scalar a, real scalar b, real scalar n)
{
	real scalar i
	real colvector out
	out = J(n, 1, .)
	if (n == 1) {
		out[1] = a
		return(out)
	}
	for (i=1; i<=n; i++) out[i] = a + (b - a) * (i - 1) / (n - 1)
	return(out)
}

real colvector lpcde_default_grid(real colvector y, real scalar ng, string scalar spacing)
{
	real scalar i
	real colvector probs, out
	if (ng < 1) ng = 19
	if (spacing == "quantile") {
		probs = lpcde_linspace(.1, .9, ng)
		out = J(ng, 1, .)
		for (i=1; i<=ng; i++) out[i] = lpcde_quantile(y, probs[i])
		return(out)
	}
	return(lpcde_linspace(lpcde_quantile(y, .1), lpcde_quantile(y, .9), ng))
}

real matrix lpcde_exponents_degree(real scalar d, real scalar degree)
{
	real scalar a
	real matrix sub, block, out
	if (d == 1) return(J(1, 1, degree))
	out = J(0, d, .)
	for (a=degree; a>=0; a--) {
		sub = lpcde_exponents_degree(d - 1, degree - a)
		block = (J(rows(sub), 1, a), sub)
		out = out \ block
	}
	return(out)
}

real matrix lpcde_exponents(real scalar d, real scalar p)
{
	real scalar deg
	real matrix out
	out = J(1, d, 0)
	for (deg=1; deg<=p; deg++) out = out \ lpcde_exponents_degree(d, deg)
	return(out)
}

real rowvector lpcde_poly_base(real rowvector x, real scalar p)
{
	real scalar i, j, deg, val
	real matrix E
	real rowvector out
	E = lpcde_exponents(cols(x), p)
	out = J(1, rows(E), 1)
	for (i=1; i<=rows(E); i++) {
		deg = sum(E[i, .])
		val = 1
		for (j=1; j<=cols(x); j++) {
			if (E[i, j] > 0) val = val * x[j]^E[i, j]
		}
		out[i] = val / lpcde_factorial(deg)
	}
	return(out)
}

real colvector lpcde_basis_vec(real rowvector x, real scalar p, real scalar mu)
{
	real matrix E
	real colvector out
	if (cols(x) == 1) {
		out = J(p + 1, 1, 0)
		out[mu + 1] = 1
		return(out)
	}
	E = lpcde_exponents(cols(x), p)
	out = (rowsum(E) :== mu)
	return(out)
}

real matrix lpcde_univar_basis(real colvector x, real scalar p)
{
	real scalar j
	real matrix R
	R = J(rows(x), p + 1, 1)
	if (p == 0) return(R)
	for (j=1; j<=p; j++) R[., j + 1] = R[., j] :* x :/ j
	return(R)
}

real scalar lpcde_kernel_scalar(real rowvector x, string scalar kernel)
{
	real scalar j, out, z
	out = 1
	for (j=1; j<=cols(x); j++) {
		z = x[j]
		if (kernel == "uniform") {
			if (abs(z) <= 1) out = out * .5
			else out = 0
		}
		else if (kernel == "triangular") {
			out = out * ((1 - abs(z)) * (abs(z) <= 1))
		}
		else {
			out = out * (.75 * (1 - z^2) * (abs(z) <= 1))
		}
	}
	return(out)
}

real colvector lpcde_kernel_vec(real colvector x, string scalar kernel)
{
	if (kernel == "uniform") return(.5 :* (abs(x) :<= 1))
	if (kernel == "triangular") return((1 :- abs(x)) :* (abs(x) :<= 1))
	return(.75 :* (1 :- x:^2) :* (abs(x) :<= 1))
}

real matrix lpcde_S_x(real matrix X, real scalar p, string scalar kernel)
{
	real scalar i, kdim, kval
	real matrix S, R
	real rowvector r
	if (rows(X) == 0) {
		kdim = rows(lpcde_basis_vec(J(1, cols(X), 0), p, 0))
		return(J(kdim, kdim, 0))
	}
	if (cols(X) == 1) {
		R = lpcde_univar_basis(X[., 1], p)
		S = J(cols(R), cols(R), 0)
		for (i=1; i<=rows(R); i++) {
			kval = lpcde_kernel_scalar(X[i, .], kernel)
			S = S + kval * (R[i, .]' * R[i, .])
		}
		return(S)
	}
	kdim = rows(lpcde_basis_vec(X[1, .], p, 0))
	S = J(kdim, kdim, 0)
	for (i=1; i<=rows(X); i++) {
		r = lpcde_poly_base(X[i, .], p)
		kval = lpcde_kernel_scalar(X[i, .], kernel)
		S = S + kval * (r' * r)
	}
	return(S)
}

real matrix lpcde_solve_checked(real matrix A, real scalar k)
{
	if (rows(A) != cols(A) | rows(A) == 0) return(J(k, k, 0))
	if (rank(A) < cols(A)) return(J(k, k, 0))
	return(luinv(A))
}

real colvector lpcde_b_x(real matrix X, real matrix S_inv, real colvector e, real scalar p, string scalar kernel)
{
	real scalar i
	real colvector out, w
	real matrix R
	real colvector K
	if (rows(X) == 0) return(J(0, 1, 0))
	w = S_inv * e
	out = J(rows(X), 1, 0)
	if (cols(X) == 1) {
		R = lpcde_univar_basis(X[., 1], p)
		K = lpcde_kernel_vec(X[., 1], kernel)
		for (i=1; i<=rows(X); i++) out[i] = (R[i, .] * K[i]) * w
		return(out)
	}
	for (i=1; i<=rows(X); i++) out[i] = (lpcde_poly_base(X[i, .], p) * lpcde_kernel_scalar(X[i, .], kernel)) * w
	return(out)
}

real matrix lpcde_fhat(real matrix x_data, real colvector y_data, real rowvector x, real colvector y_grid,
	real scalar p, real scalar q, real scalar mu, real scalar nu, real colvector h, string scalar kernel)
{
	real scalar n, d, ng, j, hval
	real matrix out, x_idx, x_sorted, x_scaled, sx_inv, sy_inv
	real colvector idx, y_idx, ord, y_sorted, y_scaled, y_elems, bx, ax
	real colvector e_nu, e_mu
	n = rows(y_data)
	d = cols(x_data)
	ng = rows(y_grid)
	out = J(ng, 2, 0)
	e_nu = lpcde_basis_vec(x, q, nu)
	e_mu = lpcde_basis_vec(J(1, 1, 0), p, mu)

	if (min(h) == max(h)) {
		hval = h[1]
		idx = lpcde_which(rowsum(abs(lpcde_sweep_row(x_data, x, "-")) :<= hval) :== d)
		if (rows(idx) == 0) return(out)
		x_idx = x_data[idx, .]
		y_idx = y_data[idx]
		ord = order(y_idx, 1)
		y_sorted = y_idx[ord]
		x_sorted = x_idx[ord, .]
		x_scaled = lpcde_sweep_row(x_sorted, x, "-") :/ (hval^d)
		sx_inv = lpcde_solve_checked(lpcde_S_x(x_scaled, q, kernel) :/ (n * hval^d), rows(e_nu))
		bx = lpcde_b_x(x_scaled, sx_inv, e_nu, q, kernel)
		for (j=1; j<=ng; j++) {
			y_scaled = (y_sorted :- y_grid[j]) :/ hval
			y_elems = lpcde_which(abs(y_scaled) :<= 1)
			out[j, 2] = rows(y_elems)
			if (rows(y_elems) <= 5) continue
			if (mu == 0) {
				sx_inv = lpcde_solve_checked(lpcde_S_x(x_scaled[y_elems, .], q, kernel) :/ (n * hval^d), rows(e_nu))
				bx = lpcde_b_x(x_scaled, sx_inv, e_nu, q, kernel)
				sy_inv = lpcde_solve_checked(lpcde_S_x(y_scaled, p, kernel) :/ (n * hval), rows(e_mu))
				ax = lpcde_b_x(y_scaled, sy_inv, e_mu, p, kernel)
				out[j, 1] = ax[y_elems]' * runningsum(bx[y_elems])
			}
			else {
				sy_inv = lpcde_solve_checked(lpcde_S_x(y_scaled[y_elems], p, kernel) :/ (n * hval), rows(e_mu))
				ax = lpcde_b_x(y_scaled[y_elems], sy_inv, e_mu, p, kernel)
				out[j, 1] = ax' * runningsum(bx[y_elems])
			}
		}
		out[., 1] = out[., 1] :/ (n^2 * hval^(d + mu + nu + 1))
		return(out)
	}

	for (j=1; j<=ng; j++) {
		hval = h[j]
		idx = lpcde_which(rowsum(abs(lpcde_sweep_row(x_data, x, "-")) :<= hval) :== d)
		if (rows(idx) == 0) continue
		x_idx = x_data[idx, .]
		y_idx = y_data[idx]
		ord = order(y_idx, 1)
		y_sorted = y_idx[ord]
		x_sorted = x_idx[ord, .]
		x_scaled = lpcde_sweep_row(x_sorted, x, "-") :/ (hval^d)
		sx_inv = lpcde_solve_checked(lpcde_S_x(x_scaled, q, kernel) :/ (n * hval^d), rows(e_nu))
		bx = lpcde_b_x(x_scaled, sx_inv, e_nu, q, kernel)
		y_scaled = (y_sorted :- y_grid[j]) :/ hval
		y_elems = lpcde_which(abs(y_scaled) :<= 1)
		out[j, 2] = rows(y_elems)
		if (rows(y_elems) <= 5) continue
		if (mu == 0) {
			sx_inv = lpcde_solve_checked(lpcde_S_x(x_scaled[y_elems, .], q, kernel) :/ (n * hval^d), rows(e_nu))
			bx = lpcde_b_x(x_scaled, sx_inv, e_nu, q, kernel)
			sy_inv = lpcde_solve_checked(lpcde_S_x(y_scaled, p, kernel) :/ (n * hval), rows(e_mu))
			ax = lpcde_b_x(y_scaled, sy_inv, e_mu, p, kernel)
			out[j, 1] = ax[y_elems]' * runningsum(bx[y_elems])
		}
		else {
			sy_inv = lpcde_solve_checked(lpcde_S_x(y_scaled[y_elems], p, kernel) :/ (n * hval), rows(e_mu))
			ax = lpcde_b_x(y_scaled[y_elems], sy_inv, e_mu, p, kernel)
			out[j, 1] = ax' * runningsum(bx[y_elems])
		}
		out[j, 1] = out[j, 1] / (n^2 * hval^(d + mu + nu + 1))
	}
	return(out)
}

real matrix lpcde_cov_hat(real matrix x_data, real colvector y_data, real rowvector x, real colvector y_grid,
	real scalar p, real scalar q, real scalar mu, real scalar nu, real colvector h, string scalar kernel, string scalar covflag)
{
	real scalar n, d, ng, i, j, k, hval, y, yp, val, kk
	real matrix c_hat, x_idx, x_sorted, x_scaled, sx_inv, sy_inv, syp_inv
	real colvector idx, y_idx, ord, y_sorted, y_scaled, yp_scaled, y_elems, yp_elems, elems
	real colvector bx, ay, ayp, aj, ak, theta, e_nu, e_mu
	n = rows(y_data)
	d = cols(x_data)
	ng = rows(y_grid)
	if (covflag == "diag") return(J(ng, ng, .))
	c_hat = J(ng, ng, 0)
	e_nu = lpcde_basis_vec(x, q, nu)
	e_mu = lpcde_basis_vec(J(1, 1, 1), p, mu)

	theta = J(ng, 1, 0)
	for (i=1; i<=ng; i++) {
		if (mu == 0) theta[i] = lpcde_fhat(x_data, y_data, x, J(1, 1, y_grid[i]), 2, 1, 0, 0, J(1, 1, h[i]), kernel)[1, 1]
		else theta[i] = lpcde_fhat(x_data, y_data, x, J(1, 1, y_grid[i]), p, q, mu, nu, J(1, 1, h[i]), kernel)[1, 1]
	}

	if (min(h) == max(h)) {
		hval = h[1]
		idx = lpcde_which(rowsum(abs(lpcde_sweep_row(x_data, x, "-")) :<= hval) :== d)
		if (rows(idx) == 0) return(c_hat)
		x_idx = x_data[idx, .]
		y_idx = y_data[idx]
		ord = order(y_idx, 1)
		y_sorted = y_idx[ord]
		x_sorted = x_idx[ord, .]
		x_scaled = lpcde_sweep_row(x_sorted, x, "-") :/ (hval^d)
		sx_inv = lpcde_solve_checked(lpcde_S_x(x_scaled, q, kernel) :/ (n * hval^d), rows(e_nu))
		bx = lpcde_b_x(x_scaled, sx_inv, e_nu, q, kernel)
		for (i=1; i<=ng; i++) {
			for (j=1; j<=i; j++) {
				y = y_grid[i]
				yp = y_grid[j]
				y_scaled = (y_sorted :- y) :/ hval
				yp_scaled = (y_sorted :- yp) :/ hval
				y_elems = lpcde_which(abs(y_scaled) :<= 1)
				yp_elems = lpcde_which(abs(yp_scaled) :<= 1)
				elems = lpcde_intersect(y_elems, yp_elems)
				val = 0
				if (rows(elems) > 5) {
					if (mu == 0) {
						sx_inv = lpcde_solve_checked(lpcde_S_x(x_scaled[y_elems, .] :/ (n * hval^d), q, kernel), rows(e_nu))
						bx = lpcde_b_x(x_scaled, sx_inv, e_nu, q, kernel)
					}
					sy_inv = lpcde_solve_checked(lpcde_S_x(y_scaled[y_elems], p, kernel) :/ (n * hval), rows(e_mu))
					syp_inv = lpcde_solve_checked(lpcde_S_x(yp_scaled[yp_elems], p, kernel) :/ (n * hval), rows(e_mu))
					ay = lpcde_b_x(y_scaled[y_elems], sy_inv, e_mu, p, kernel)
					ayp = lpcde_b_x(yp_scaled[elems], syp_inv, e_mu, p, kernel)
					aj = runningsum(ay)
					ak = runningsum(ayp)
					for (k=1; k<=rows(elems); k++) {
						kk = k
						val = val + bx[elems[k]]^2 * (aj[k]*ak[k] + (n-kk)*ayp[k]*aj[k] + (n-kk)*ay[k]*ak[k] + (n-kk)^2*ay[k]*ayp[k])
					}
				}
				val = val / (n * (n - 1)^2) - theta[i] * theta[j] / n^2
				c_hat[i, j] = val
				c_hat[j, i] = val
			}
		}
		if (mu == 0) c_hat = c_hat :* (1 / (n * hval^(d + mu + nu)))^2
		else {
			c_hat = c_hat :* (1 / (n * hval^(d + mu + nu + 1)))^2
			for (i=1; i<=ng; i++) c_hat[i, i] = c_hat[i, i] / 2
		}
		return(c_hat)
	}

	for (i=1; i<=ng; i++) {
		for (j=1; j<=i; j++) {
			hval = (h[i] + h[j]) / 2
			idx = lpcde_which(rowsum(abs(lpcde_sweep_row(x_data, x, "-")) :<= hval) :== d)
			if (rows(idx) == 0) continue
			x_idx = x_data[idx, .]
			y_idx = y_data[idx]
			ord = order(y_idx, 1)
			y_sorted = y_idx[ord]
			x_sorted = x_idx[ord, .]
			x_scaled = lpcde_sweep_row(x_sorted, x, "-") :/ (hval^d)
			sx_inv = lpcde_solve_checked(lpcde_S_x(x_scaled, q, kernel) :/ (n * hval^d), rows(e_nu))
			bx = lpcde_b_x(x_scaled, sx_inv, e_nu, q, kernel)
			y_scaled = (y_sorted :- y_grid[i]) :/ h[i]
			yp_scaled = (y_sorted :- y_grid[j]) :/ h[j]
			y_elems = lpcde_which(abs(y_scaled) :<= 1)
			yp_elems = lpcde_which(abs(yp_scaled) :<= 1)
			elems = lpcde_intersect(y_elems, yp_elems)
			val = 0
			if (rows(elems) > 5) {
				if (mu == 0) {
					sx_inv = lpcde_solve_checked(lpcde_S_x(x_scaled[y_elems, .] :/ (n * hval^d), q, kernel), rows(e_nu))
					bx = lpcde_b_x(x_scaled, sx_inv, e_nu, q, kernel)
				}
				sy_inv = lpcde_solve_checked(lpcde_S_x(y_scaled[y_elems], p, kernel) :/ (n * h[i]), rows(e_mu))
				syp_inv = lpcde_solve_checked(lpcde_S_x(yp_scaled[yp_elems], p, kernel) :/ (n * h[j]), rows(e_mu))
				ay = lpcde_b_x(y_scaled[elems], sy_inv, e_mu, p, kernel)
				ayp = lpcde_b_x(yp_scaled[elems], syp_inv, e_mu, p, kernel)
				aj = runningsum(ay)
				ak = runningsum(ayp)
				for (k=1; k<=rows(elems); k++) {
					kk = k
					val = val + bx[elems[k]]^2 * (aj[k]*ak[k] + (n-kk)*ayp[k]*aj[k] + (n-kk)*ay[k]*ak[k] + (n-kk)^2*ay[k]*ayp[k])
				}
			}
			val = val / (n * (n - 1)^2) - theta[i] * theta[j] / n^2
			c_hat[i, j] = val
			c_hat[j, i] = val
		}
	}
	for (i=1; i<=ng; i++) {
		for (j=1; j<=ng; j++) {
			if (mu == 0) c_hat[i, j] = c_hat[i, j] / (n * h[i]^(d + mu + nu)) / (n * h[j]^(d + mu + nu))
			else c_hat[i, j] = c_hat[i, j] / (n * h[i]^(d + mu + nu + 1)) / (n * h[j]^(d + mu + nu + 1))
		}
		if (mu != 0) c_hat[i, i] = c_hat[i, i] / 2
	}
	return(c_hat)
}

real scalar lpcde_int_val(real scalar l, real scalar a, real scalar b, string scalar kernel)
{
	real scalar num, num1, num2, denom1, denom2
	if (kernel == "triangular") {
		if (a >= 0 & b >= 0) {
			num = a^(l+1) * (-2 + a + (-1+a)*l) + b^(l+1) * (2+l-b*(1+l))
			return(num / ((l+1)*(l+2)))
		}
		if (a < 0 & b > 0) {
			num1 = -a^(l+1) * (2+a+l+a*l)
			denom1 = 2 + 3*l + l^2
			num2 = (2+l)*b^(l+1) - (l+1)*b^(l+2)
			denom2 = (l+1)*(l+2)
			return(num1/denom1 + num2/denom2)
		}
		num = -a^(l+1) * (2+a+l+a*l) + b^(l+1) * (2+l+b+b*l)
		return(num / ((l+1)*(l+2)))
	}
	if (kernel == "uniform") return(.5 * (b^(l+1) - a^(l+1)) / (l+1))
	num = (l+1)*a^(l+3) - (l+3)*a^(l+1) + b^(l+1)*(3+l-b^2*(l+1))
	return(.75 * num / ((l+1)*(l+3)))
}

real matrix lpcde_S_exact(real scalar lower, real scalar upper, real scalar eval_pt, real scalar p, string scalar kernel)
{
	real scalar i, j, a, b
	real colvector poly
	real matrix S
	a = max((lower, -1))
	b = min((upper, 1))
	poly = J(2*p + 1, 1, 0)
	if (a < b) {
		for (i=0; i<=2*p; i++) poly[i+1] = lpcde_int_val(i, a, b, kernel)
	}
	S = J(p+1, p+1, 0)
	for (j=1; j<=p+1; j++) S[., j] = poly[j..(j+p)]
	for (i=0; i<=p; i++) {
		for (j=0; j<=p; j++) S[i+1, j+1] = normalden(eval_pt) * S[i+1, j+1] / (lpcde_factorial(i) * lpcde_factorial(j))
	}
	return(S)
}

real colvector lpcde_c_exact(real scalar lower, real scalar upper, real scalar eval_pt, real scalar m, real scalar p, string scalar kernel)
{
	real scalar i, a, b
	real colvector v
	a = max((lower, -1))
	b = min((upper, 1))
	v = J(p+1, 1, 0)
	if (a < b) {
		for (i=0; i<=p; i++) v[i+1] = lpcde_int_val(i + m, a, b, kernel) / lpcde_factorial(m)
	}
	return(v)
}

real matrix lpcde_T_y_exact(real scalar lower, real scalar upper, real scalar p)
{
	real scalar i, j, a, b, num1, num2, num3, denom1, denom2, denom3
	real matrix T
	a = max((lower, -1))
	b = min((upper, 1))
	T = J(p+1, p+1, 0)
	if (a < b) {
		for (i=0; i<=p; i++) {
			for (j=0; j<=p; j++) {
				num1 = -(b^(i+j+3) - a^(i+j+3))
				denom1 = (i+1)*(i+2)*(i+j+3)
				num2 = -(a^(i+2)*(b^(j+1)-a^(j+1)))
				denom2 = (i+2)*(j+1)
				num3 = b^(i+1)*(b^(j+2)-a^(j+2))
				denom3 = (i+2)*(j+2)
				T[i+1, j+1] = (num1/denom1 + num2/denom2 + num3/denom3) / 4 / (lpcde_factorial(i) * lpcde_factorial(j))
			}
		}
	}
	return(T)
}

real matrix lpcde_T_x(real matrix x_data, real scalar x, real scalar q, real scalar h, string scalar kernel)
{
	real matrix R, T
	real colvector z, K
	real scalar i
	z = (x_data[., 1] :- x) :/ h
	R = lpcde_univar_basis(z, q)
	K = lpcde_kernel_vec(z, kernel)
	T = J(cols(R), cols(R), 0)
	for (i=1; i<=rows(R); i++) T = T + K[i]^2 * (R[i, .]' * R[i, .])
	return(T / rows(x_data))
}

real scalar lpcde_hermite(real scalar n, real scalar z)
{
	real scalar h0, h1, h2, i
	if (n == 0) return(1)
	if (n == 1) return(z)
	h0 = 1
	h1 = z
	for (i=2; i<=n; i++) {
		h2 = z*h1 - (i-1)*h0
		h0 = h1
		h1 = h2
	}
	return(h1)
}

real scalar lpcde_normal_dgps(real scalar x, real scalar v, real scalar mu, real scalar sd)
{
	real scalar z, order
	if (v == 0) return(normal((x - mu) / sd))
	z = (x - mu) / sd
	order = v - 1
	return(((-1)^order) * lpcde_hermite(order, z) * normalden(z) / (sd^v))
}

real colvector lpcde_bw_d1(real colvector y_data, real matrix x_data, real colvector y_grid, real rowvector x,
	real scalar p, real scalar q, real scalar mu, real scalar nu, string scalar kernel, real scalar regularize, real scalar integrated)
{
	real scalar n, ng, sd_y, sd_x, mx, my, bx, lower_x, upper_x, alpha, hscalar, id
	real scalar j, y, lower_y, upper_y, cdf_hat
	real matrix bias, Sx, Sy, Tx, Ty
	real colvector v, h, cx, cy
	n = rows(y_data)
	ng = rows(y_grid)
	sd_y = lpcde_sd(y_data)
	sd_x = lpcde_sd(x_data[., 1])
	mx = mean(x_data[., 1])
	my = mean(y_data)
	if (integrated) bx = .5
	else bx = 1.06 * n^(-1/5) * sd_x
	lower_x = min(x_data[., 1]) - x[1]
	upper_x = max(x_data[., 1]) - x[1]
	bias = J(ng, 3, .)
	for (j=1; j<=ng; j++) {
		y = y_grid[j]
		lower_y = min(y_data) - y
		upper_y = max(y_data) - y
		bias[j, 1] = lpcde_normal_dgps(y, mu, my, sd_y) * lpcde_normal_dgps(x[1], 2, mx, sd_x)
		bias[j, 2] = lpcde_normal_dgps(y, p + 1, my, sd_y) * lpcde_normal_dgps(x[1], 0, mx, sd_x)
		Sx = luinv(lpcde_S_exact(lower_x, upper_x, x[1], q, kernel))
		Sy = luinv(lpcde_S_exact(lower_y, upper_y, y, p, kernel))
		cx = lpcde_c_exact(lower_x, upper_x, x[1], q + 1, q, kernel)
		cy = lpcde_c_exact(lower_y, upper_y, y, p + 1, p, kernel)
		bias[j, 1] = bias[j, 1] * (cx' * Sx)[1, nu + 1]
		bias[j, 2] = bias[j, 2] * (cy' * Sy)[1, mu + 1]
		bias[j, 3] = (bias[j, 1] + bias[j, 2])^2
	}
	v = J(ng, 1, .)
	for (j=1; j<=ng; j++) {
		y = y_grid[j]
		lower_y = min(y_data) - y
		upper_y = max(y_data) - y
		Sx = luinv(lpcde_S_exact(lower_x, upper_x, x[1], q, kernel))
		if (mu > 0) {
			Sy = luinv(lpcde_S_exact(lower_y, upper_y, y, p, kernel))
			Ty = lpcde_T_y_exact(lower_y, upper_y, p)
			if (mu == 1) Tx = lpcde_T_x(x_data, x[1], q, bx, kernel) :/ bx
			else Tx = lpcde_T_y_exact(lower_x, upper_x, q)
			v[j] = normalden(y) * normal(x[1]) * (Sy * Ty * Sy)[mu + 1, mu + 1] * (Sx * Tx * Sx)[nu + 1, nu + 1]
		}
		else {
			Tx = lpcde_T_x(x_data, x[1], q, bx, kernel) :/ bx
			cdf_hat = normal(y) * normal(x[1])
			v[j] = cdf_hat * (1 - cdf_hat) * (Sx * Tx * Sx)[nu + 1, nu + 1]
		}
	}
	if (mu == 0) alpha = 1 + 2*min((p, q)) + 2*min((mu, nu)) + 2*nu + 2
	else alpha = 1 + 2*min((p, q)) + 2*max((mu, nu)) + 1
	if (integrated) {
		hscalar = (abs(mean(v) / (2 * mean(bias[., 3]))))^(1/alpha) * n^(-1/alpha)
		hscalar = sd_y * sd_x * hscalar
		if (regularize) {
			id = min((n, 20 + q + 1))
			hscalar = max((hscalar, sort(abs(x_data[., 1] :- x[1]), 1)[id]))
			id = min((n, 20 + p + 1))
			for (j=1; j<=ng; j++) hscalar = max((hscalar, sort(abs(y_data :- y_grid[j]), 1)[id]))
		}
		return(J(ng, 1, hscalar))
	}
	h = (abs(v :/ bias[., 3])):^(1/alpha) :* n^(-1/alpha)
	if (regularize) {
		for (j=1; j<=ng; j++) {
			id = min((n, 20 + q + 4))
			h[j] = max((h[j], sort(abs(x_data[., 1] :- x[1]), 1)[id]))
			id = min((n, 20 + p + 4))
			h[j] = max((h[j], sort(abs(y_data :- y_grid[j]), 1)[id]))
		}
	}
	return(h)
}

real colvector lpcde_bw_select(real colvector y_raw, real matrix x_raw, real rowvector x_raw_eval, real colvector y_grid,
	real scalar p, real scalar q, real scalar mu, real scalar nu, string scalar kernel, string scalar bwselect, real scalar regularize)
{
	real scalar n, d, j, h0
	real rowvector sd_x, mx, x_scaled
	real colvector y_scaled, h
	real matrix x_scaled_data
	n = rows(y_raw)
	d = cols(x_raw)
	sd_x = lpcde_colsd(x_raw)
	mx = lpcde_colmean(x_raw)
	y_scaled = y_raw :/ lpcde_sd(y_raw)
	x_scaled_data = lpcde_sweep_row(x_raw, sd_x, "/")
	x_scaled = (x_raw_eval - mx) :/ sd_x
	if (d == 1) return(lpcde_bw_d1(y_scaled, x_scaled_data, y_grid, x_scaled, p, q, mu, nu, kernel, regularize, bwselect == "imse-rot"))

	h0 = lpcde_sd(y_raw) * lpcde_rowmean(sd_x) * n^(-1 / (d + 2*min((p, q)) + 2*max((mu, nu)) + 1))
	h = J(rows(y_grid), 1, h0)
	if (regularize) {
		for (j=1; j<=rows(y_grid); j++) {
			h[j] = max((h[j], sort(abs(y_raw :- y_grid[j]), 1)[min((n, 20 + p + 1))]))
		}
	}
	return(h)
}

real matrix lpcde_bw_result(real colvector y_raw, real matrix x_raw, real rowvector x_raw_eval, real colvector y_grid,
	real scalar p, real scalar q, real scalar mu, real scalar nu, string scalar kernel, string scalar bwselect, real scalar regularize)
{
	real scalar i, d
	real rowvector sd_x, mx, x_scaled
	real colvector y_scaled, bw, idx
	real matrix x_scaled_data, out
	d = cols(x_raw)
	sd_x = lpcde_colsd(x_raw)
	mx = lpcde_colmean(x_raw)
	y_scaled = y_raw :/ lpcde_sd(y_raw)
	x_scaled_data = lpcde_sweep_row(x_raw, sd_x, "/")
	x_scaled = (x_raw_eval - mx) :/ sd_x
	bw = lpcde_bw_select(y_raw, x_raw, x_raw_eval, y_grid, p, q, mu, nu, kernel, bwselect, regularize)
	out = J(rows(y_grid), 4, .)
	out[., 1] = y_grid
	out[., 2] = bw
	for (i=1; i<=rows(y_grid); i++) {
		idx = lpcde_which(rowsum(abs(lpcde_sweep_row(x_scaled_data, x_scaled, "-")) :<= bw[i]) :== d)
		if (rows(idx) == 0) out[i, 3] = 0
		else out[i, 3] = sum(abs(y_scaled[idx] :- y_grid[i]) :<= bw[i])
		out[i, 4] = i
	}
	return(out)
}

real matrix lpcde_estimate(real colvector y_raw, real matrix x_raw, real rowvector x_raw_eval, real colvector y_grid,
	real colvector bw, real scalar p, real scalar q, real scalar prbc, real scalar qrbc, real scalar mu, real scalar nu,
	real scalar rbc, string scalar kernel, string scalar covflag, real scalar nonneg, real scalar normalize, real scalar level)
{
	real scalar ng, d, sd_y, i, z, total
	real rowvector sd_x, mx, x_scaled
	real matrix x_scaled_data, fit, fit_rbc, cov, cov_rbc, out
	real colvector se, se_rbc, grid_diff
	ng = rows(y_grid)
	d = cols(x_raw)
	sd_y = lpcde_sd(y_raw)
	sd_x = lpcde_colsd(x_raw)
	mx = lpcde_colmean(x_raw)
	x_scaled_data = lpcde_sweep_row(lpcde_sweep_row(x_raw, mx, "-"), sd_x, "/")
	x_scaled = (x_raw_eval - mx) :/ sd_x
	fit = lpcde_fhat(x_scaled_data, y_raw, x_scaled, y_grid, p, q, mu, nu, bw, kernel)
	if (covflag == "off") {
		cov = J(ng, ng, .)
		se = J(ng, 1, .)
	}
	else {
		cov = lpcde_cov_hat(x_scaled_data, y_raw, x_scaled, y_grid, p, q, mu, nu, bw, kernel, covflag)
		se = sqrt(abs(diagonal(cov))) :* sd_y :* lpcde_rowmean(sd_x)
	}
	if (rbc & (prbc != p | qrbc != q)) {
		fit_rbc = lpcde_fhat(x_scaled_data, y_raw, x_scaled, y_grid, prbc, qrbc, mu, nu, bw, kernel)
		if (covflag == "off") {
			cov_rbc = J(ng, ng, .)
			se_rbc = J(ng, 1, .)
		}
		else {
			cov_rbc = lpcde_cov_hat(x_scaled_data, y_raw, x_scaled, y_grid, prbc, qrbc, mu, nu, bw, kernel, covflag)
			se_rbc = sqrt(abs(diagonal(cov_rbc))) :* sd_y :* lpcde_rowmean(sd_x)
		}
	}
	else {
		fit_rbc = fit
		cov_rbc = cov
		se_rbc = se
	}
	if (nonneg) {
		for (i=1; i<=ng; i++) if (fit[i, 1] < 0) fit[i, 1] = 0
	}
	if (normalize & ng > 1) {
		grid_diff = (y_grid[2..ng] - y_grid[1..(ng-1)]) \ (y_grid[ng] - y_grid[ng-1])
		total = sum(fit[., 1] :* grid_diff)
		if (total != 0) fit[., 1] = fit[., 1] :/ total
	}
	z = invnormal(1 - (1 - level/100) / 2)
	out = J(ng, 10, .)
	out[., 1] = y_grid
	out[., 2] = bw
	out[., 3] = fit[., 2]
	out[., 4] = fit[., 1]
	out[., 5] = fit_rbc[., 1]
	out[., 6] = se
	out[., 7] = se_rbc
	if (rbc) {
		out[., 8] = fit_rbc[., 1] :- z :* se_rbc
		out[., 9] = fit_rbc[., 1] :+ z :* se_rbc
	}
	else {
		out[., 8] = fit[., 1] :- z :* se
		out[., 9] = fit[., 1] :+ z :* se
	}
	out[., 10] = (1::ng)
	st_matrix("CovMat", cov)
	st_matrix("CovMat_RBC", cov_rbc)
	return(out)
}

real matrix lpcde_select_result(real matrix out, string scalar rgrid_name, string scalar rindex_name)
{
	real scalar i
	real colvector rvals, sel
	if (rgrid_name != "") {
		rvals = lpcde_nonmissing(st_data(., rgrid_name))
		sel = J(rows(rvals), 1, .)
		for (i=1; i<=rows(rvals); i++) sel[i] = order(abs(out[., 1] :- rvals[i]), 1)[1]
		return(out[sel, .])
	}
	if (rindex_name != "") {
		rvals = lpcde_nonmissing(st_data(., rindex_name))
		sel = J(rows(rvals), 1, .)
		for (i=1; i<=rows(rvals); i++) sel[i] = order(abs(out[., cols(out)] :- rvals[i]), 1)[1]
		return(out[sel, .])
	}
	return(out)
}

void lpcde_stata_run(string scalar varlist, string scalar touse, string scalar grid_name, string scalar rgrid_name, string scalar rindex_name,
	string scalar bw_name, real scalar has_bw, real scalar bw_is_var, real scalar bw_scalar, real scalar has_xeval,
	real scalar ng, string scalar gridspacing, real scalar p, real scalar q, real scalar prbc, real scalar qrbc,
	real scalar mu, real scalar nu, real scalar rbc, string scalar kernel, string scalar bwselect, string scalar covflag,
	real scalar nonneg, real scalar normalize, real scalar level)
{
	real scalar j
	real matrix data, x_raw, out
	real colvector y_raw, y_grid, bw
	real rowvector x_eval
	data = st_data(., varlist, touse)
	y_raw = data[., 1]
	x_raw = data[., 2..cols(data)]
	if (has_xeval) x_eval = st_matrix("Xeval")
	else {
		x_eval = J(1, cols(x_raw), .)
		for (j=1; j<=cols(x_raw); j++) x_eval[j] = lpcde_quantile(x_raw[., j], .5)
	}
	if (grid_name != "") y_grid = lpcde_nonmissing(st_data(., grid_name))
	else y_grid = lpcde_default_grid(y_raw, ng, gridspacing)
	if (has_bw) {
		if (bw_is_var) bw = lpcde_nonmissing(st_data(., bw_name))
		else bw = J(rows(y_grid), 1, bw_scalar)
	}
	else bw = lpcde_bw_select(y_raw, x_raw, x_eval, y_grid, p, q, mu, nu, kernel, bwselect, 1)
	if (rows(bw) == 1) bw = J(rows(y_grid), 1, bw[1])
	out = lpcde_estimate(y_raw, x_raw, x_eval, y_grid, bw, p, q, prbc, qrbc, mu, nu, rbc, kernel, covflag, nonneg, normalize, level)
	out = lpcde_select_result(out, rgrid_name, rindex_name)
	st_matrix("Result", out)
}

void lpcde_stata_bw_run(string scalar varlist, string scalar touse, string scalar grid_name, string scalar rgrid_name, string scalar rindex_name,
	real scalar has_xeval, real scalar ng, string scalar gridspacing, real scalar p, real scalar q, real scalar mu, real scalar nu,
	string scalar kernel, string scalar bwselect, real scalar regularize)
{
	real scalar j
	real matrix data, x_raw, out
	real colvector y_raw, y_grid
	real rowvector x_eval
	data = st_data(., varlist, touse)
	y_raw = data[., 1]
	x_raw = data[., 2..cols(data)]
	if (has_xeval) x_eval = st_matrix("Xeval")
	else {
		x_eval = J(1, cols(x_raw), .)
		for (j=1; j<=cols(x_raw); j++) x_eval[j] = lpcde_quantile(x_raw[., j], .5)
	}
	if (grid_name != "") y_grid = lpcde_nonmissing(st_data(., grid_name))
	else y_grid = lpcde_default_grid(y_raw :/ lpcde_sd(y_raw), ng, gridspacing)
	out = lpcde_bw_result(y_raw, x_raw, x_eval, y_grid, p, q, mu, nu, kernel, bwselect, regularize)
	out = lpcde_select_result(out, rgrid_name, rindex_name)
	st_matrix("Result", out)
}

end
