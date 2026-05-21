from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
from typing import Any
import warnings

import numpy as np
import pandas as pd
from scipy.special import eval_hermitenorm
from scipy.stats import multivariate_normal, norm


ArrayLike = Any


@dataclass
class LPCDEResult:
    estimate: pd.DataFrame
    cov_mat: dict[str, np.ndarray | float]
    eff_n: np.ndarray
    opt: dict[str, Any]
    singular_flag: bool = False

    @property
    def Estimate(self) -> pd.DataFrame:
        return self.estimate

    @property
    def CovMat(self) -> dict[str, np.ndarray | float]:
        return self.cov_mat

    def coef(self, rbc: bool = False) -> pd.DataFrame:
        col = "est_RBC" if rbc else "est"
        return self.estimate[["y_grid", "bw", col]].copy()

    def vcov(self, rbc: bool = False) -> np.ndarray | float:
        key = "CovMat_RBC" if rbc else "CovMat"
        return self.cov_mat[key]

    def confint(self, level: float = 0.95, rbc: bool = True) -> pd.DataFrame:
        alpha = 1.0 - level
        z = norm.ppf(1.0 - alpha / 2.0)
        center = "est_RBC" if rbc else "est"
        se = "se_RBC" if rbc else "se"
        out = self.estimate[["y_grid", center, se]].copy()
        out["ci_l"] = out[center] - z * out[se]
        out["ci_r"] = out[center] + z * out[se]
        return out

    def summary(self) -> pd.DataFrame:
        return self.estimate.copy()

    def plot(self, rbc: bool = True, ci: bool = True, ax: Any | None = None) -> Any:
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError("matplotlib is required for LPCDEResult.plot().") from exc

        if ax is None:
            _, ax = plt.subplots()
        center = "est_RBC" if rbc else "est"
        se = "se_RBC" if rbc else "se"
        x = self.estimate["y_grid"].to_numpy(dtype=float)
        y = self.estimate[center].to_numpy(dtype=float)
        ax.plot(x, y, label=center)
        if ci and not np.isnan(self.estimate[se].to_numpy(dtype=float)).all():
            z = norm.ppf(0.975)
            err = z * self.estimate[se].to_numpy(dtype=float)
            ax.fill_between(x, y - err, y + err, alpha=0.2, label="95% CI")
        ax.set_xlabel("y_grid")
        ax.set_ylabel("estimate")
        ax.legend()
        return ax


@dataclass
class LPBWCDEResult:
    bw: pd.DataFrame
    opt: dict[str, Any]

    @property
    def BW(self) -> pd.DataFrame:
        return self.bw

    def coef(self) -> pd.DataFrame:
        return self.bw.copy()

    def summary(self) -> pd.DataFrame:
        return self.bw.copy()


def coef(obj: LPCDEResult | LPBWCDEResult, *args: Any, **kwargs: Any) -> pd.DataFrame:
    return obj.coef(*args, **kwargs)


def vcov(obj: LPCDEResult, *args: Any, **kwargs: Any) -> np.ndarray | float:
    return obj.vcov(*args, **kwargs)


def confint(obj: LPCDEResult, *args: Any, **kwargs: Any) -> pd.DataFrame:
    return obj.confint(*args, **kwargs)


def summary(obj: LPCDEResult | LPBWCDEResult) -> pd.DataFrame:
    return obj.summary()


def plot(obj: LPCDEResult, *args: Any, **kwargs: Any) -> Any:
    return obj.plot(*args, **kwargs)


def lpcde(
    x_data: ArrayLike,
    y_data: ArrayLike,
    y_grid: ArrayLike | None = None,
    x: ArrayLike | None = None,
    bw: ArrayLike | None = None,
    p: int | None = None,
    q: int | None = None,
    p_RBC: int | None = None,
    q_RBC: int | None = None,
    mu: int | None = None,
    nu: int | None = None,
    rbc: bool = True,
    ng: int | None = None,
    cov_flag: str = "full",
    normalize: bool = False,
    nonneg: bool = False,
    grid_spacing: str = "",
    kernel_type: str = "epanechnikov",
    bw_type: str | None = None,
) -> LPCDEResult:
    y_data = _as_vector(y_data)
    x_data = _as_matrix(x_data)
    if len(y_data) != x_data.shape[0]:
        raise ValueError("x_data and y_data must have the same number of observations.")
    if np.isnan(y_data).any() or np.isnan(x_data).any():
        keep = ~np.isnan(y_data) & ~np.isnan(x_data).any(axis=1)
        warnings.warn(f"{np.sum(~keep)} missing observations are ignored.", stacklevel=2)
        y_data = y_data[keep]
        x_data = x_data[keep, :]
    if y_data.size == 0 or x_data.size == 0:
        raise ValueError("Data should be numeric, and cannot be empty.")

    cov_flag = (cov_flag or "full").lower()
    if cov_flag not in {"full", "diag", "off"}:
        raise ValueError("Incorrect covariance estimation flag provided.")

    if y_grid is None:
        if grid_spacing == "quantile":
            if bw is not None and np.size(bw) >= 2:
                y_grid = _r_quantile(y_data, np.linspace(0.1, 0.9, np.size(bw)))
            elif ng is not None:
                y_grid = _r_quantile(y_data, np.linspace(0.1, 0.9, ng))
            else:
                y_grid = _r_quantile(y_data, np.linspace(0.1, 0.9, 19))
        else:
            gmin, gmax = _r_quantile(y_data, [0.1, 0.9])
            if bw is not None and np.size(bw) >= 2:
                y_grid = np.linspace(gmin, gmax, np.size(bw))
            elif ng is not None:
                y_grid = np.linspace(gmin, gmax, ng)
            else:
                y_grid = np.linspace(gmin, gmax, 19)
    y_grid = _as_vector(y_grid)
    ng = len(y_grid)

    if x is None:
        x = np.median(x_data, axis=0)
    x = _as_vector(x)

    if p is None:
        p = 2 if mu is None else mu + 1
    if q is None:
        q = 1 if nu is None else nu + 1
    if p_RBC is None:
        p_RBC = p + 1
    if q_RBC is None:
        q_RBC = q + 1
    if mu is None:
        mu = min(1, p)
    if nu is None:
        nu = min(0, q)
    _check_order("p", p)
    _check_order("q", q)
    _check_order("p_RBC", p_RBC)
    _check_order("q_RBC", q_RBC)
    _check_order("mu", mu)
    _check_order("nu", nu)
    if p_RBC < p or q_RBC < q:
        raise ValueError("RBC polynomial order must be at least the estimation order.")
    if mu > p or nu > q:
        raise ValueError("Derivative order cannot exceed polynomial order.")

    kernel_type = _check_kernel(kernel_type)

    if bw is None:
        selected_type = "imse-rot" if bw_type is None else bw_type.lower()
        bw_result = lpbwcde(
            y_data=y_data,
            x_data=x_data,
            x=x,
            y_grid=y_grid,
            p=p,
            q=q,
            mu=mu,
            nu=nu,
            kernel_type=kernel_type,
            bw_type=selected_type,
        )
        bw_vec = bw_result.bw["bw"].to_numpy(dtype=float)
        bw_type_out = selected_type
    else:
        bw_vec = _as_vector(bw)
        if bw_vec.size == 1:
            if bw_vec[0] <= 0:
                raise ValueError("Bandwidth incorrectly specified.")
            bw_vec = np.repeat(bw_vec[0], ng)
        elif bw_vec.size != ng:
            raise ValueError("Bandwidth has to be the same length as grid.")
        if np.any(bw_vec <= 0):
            raise ValueError("Bandwidth incorrectly specified.")
        bw_type_out = "user provided"

    lpcdest = _lpcde_fn(
        y_data=y_data,
        x_data=x_data,
        y_grid=y_grid,
        x=x,
        p=p,
        q=q,
        p_RBC=p_RBC,
        q_RBC=q_RBC,
        bw=bw_vec,
        mu=mu,
        nu=nu,
        cov_flag=cov_flag,
        kernel_type=kernel_type,
        rbc=rbc,
    )
    estimate = lpcdest["est"].copy()
    if nonneg:
        estimate.loc[estimate["est"] < 0, "est"] = 0.0
    if normalize:
        grid_diff = np.r_[np.diff(y_grid), y_grid[-1] - y_grid[-2]]
        total = float(np.sum(estimate["est"].to_numpy() * grid_diff))
        estimate["est"] = estimate["est"] / total

    result = LPCDEResult(
        estimate=estimate,
        cov_mat=lpcdest["CovMat"],
        eff_n=lpcdest["eff_n"],
        singular_flag=bool(lpcdest["singular_flag"]),
        opt={
            "p": p,
            "q": q,
            "p_RBC": p_RBC,
            "q_RBC": q_RBC,
            "mu": mu,
            "nu": nu,
            "kernel": kernel_type,
            "n": len(y_data),
            "ng": ng,
            "bw_type": bw_type_out,
            "bw": bw_vec,
            "xeval": x,
            "cov_flag": cov_flag,
            "y_data_min": float(np.min(y_data)),
            "y_data_max": float(np.max(y_data)),
            "x_data_min": float(np.min(x_data)),
            "x_data_max": float(np.max(x_data)),
            "grid_min": float(np.min(y_grid)),
            "grid_max": float(np.max(y_grid)),
        },
    )
    if np.any(result.eff_n <= 5):
        warnings.warn("Some evaluation points do not have enough data to produce reliable results.", stacklevel=2)
    if result.singular_flag:
        warnings.warn("Singular matrices encountered. May affect estimates.", stacklevel=2)
    return result


def lpbwcde(
    y_data: ArrayLike,
    x_data: ArrayLike,
    x: ArrayLike | None = None,
    y_grid: ArrayLike | None = None,
    p: int | None = None,
    q: int | None = None,
    grid_spacing: str = "",
    ng: int | None = None,
    mu: int | None = None,
    nu: int | None = None,
    kernel_type: str = "epanechnikov",
    bw_type: str = "imse-rot",
    regularize: bool | None = None,
) -> LPBWCDEResult:
    y_data = _as_vector(y_data)
    x_data = _as_matrix(x_data)
    if len(y_data) != x_data.shape[0]:
        raise ValueError("x_data and y_data must have the same number of observations.")
    if np.isnan(y_data).any() or np.isnan(x_data).any():
        keep = ~np.isnan(y_data) & ~np.isnan(x_data).any(axis=1)
        warnings.warn(f"{np.sum(~keep)} missing observations are ignored.", stacklevel=2)
        y_data = y_data[keep]
        x_data = x_data[keep, :]

    d = x_data.shape[1]
    n = len(y_data)
    sd_y = _sd(y_data)
    sd_x = np.apply_along_axis(_sd, 0, x_data)
    mx = np.mean(x_data, axis=0)
    y_data_scaled = y_data / sd_y
    x_data_scaled = _r_recycle_binary(x_data, sd_x, op="/")
    if x is None:
        x_scaled = np.median(x_data_scaled, axis=0)
    else:
        x_scaled = (_as_vector(x) - mx) / sd_x

    if y_grid is None:
        if grid_spacing == "quantile":
            y_grid = _r_quantile(y_data_scaled, np.linspace(0.1, 0.9, ng or 19))
        else:
            gmin, gmax = _r_quantile(y_data_scaled, [0.1, 0.9])
            y_grid = np.linspace(gmin, gmax, ng or 19)
    y_grid = _as_vector(y_grid)
    ng = len(y_grid)

    if p is None:
        p = 2 if mu is None else mu + 1
    if q is None:
        q = 1 if nu is None else nu + 1
    if mu is None:
        mu = min(1, p)
    if nu is None:
        nu = min(0, q)
    _check_order("p", p)
    _check_order("q", q)
    _check_order("mu", mu)
    _check_order("nu", nu)
    if mu > p or nu > q:
        raise ValueError("Derivative order cannot exceed polynomial order.")

    bw_type = (bw_type or "imse-rot").lower()
    if bw_type not in {"mse-rot", "imse-rot"}:
        raise ValueError("Incorrect bandwidth selection method specified.")
    kernel_type = _check_kernel(kernel_type)
    regularize = True if regularize is None else bool(regularize)

    if bw_type == "mse-rot":
        bw = _bw_rot(y_data_scaled, x_data_scaled, y_grid, x_scaled, p, q, mu, nu, kernel_type, regularize)
    else:
        bw = _bw_irot(y_data_scaled, x_data_scaled, y_grid, x_scaled, p, q, mu, nu, kernel_type, regularize)
    bw = np.asarray(bw, dtype=float).reshape(-1)
    if bw.size == 1:
        bw = np.repeat(bw[0], ng)

    bw_frame = pd.DataFrame({"y_grid": y_grid, "bw": bw, "nh": np.nan})
    for i in range(ng):
        if d == 1:
            x_idx = np.flatnonzero(np.abs(x_data_scaled[:, 0] - x_scaled[0]) <= bw[i])
        else:
            x_idx = np.flatnonzero(np.sum(np.abs(x_data_scaled - x_scaled) <= bw[i], axis=1) == d)
        bw_frame.loc[i, "nh"] = np.sum(np.abs(y_data_scaled[x_idx] - y_grid[i]) <= bw[i])

    return LPBWCDEResult(
        bw=bw_frame,
        opt={
            "x": x_scaled * sd_x + mx,
            "p": p,
            "q": q,
            "mu": mu,
            "nu": nu,
            "kernel_type": kernel_type,
            "n": n,
            "ng": ng,
            "bw_type": bw_type,
            "data_min": float(np.min(y_data_scaled)),
            "data_max": float(np.max(y_data_scaled)),
            "grid_min": float(np.min(y_grid)),
            "grid_max": float(np.max(y_grid)),
        },
    )


def _lpcde_fn(
    y_data: np.ndarray,
    x_data: np.ndarray,
    y_grid: np.ndarray,
    x: np.ndarray,
    p: int,
    q: int,
    p_RBC: int,
    q_RBC: int,
    bw: np.ndarray,
    mu: int,
    nu: int,
    cov_flag: str,
    kernel_type: str,
    rbc: bool,
) -> dict[str, Any]:
    sd_y = _sd(y_data)
    sd_x = np.apply_along_axis(_sd, 0, x_data)
    mx = np.mean(x_data, axis=0)
    d = x_data.shape[1]
    x_data_scaled = _r_recycle_binary(x_data - mx.reshape(1, d), sd_x, op="/")
    x_scaled = (_as_vector(x).reshape(1, d) - mx.reshape(1, d)) / sd_x.reshape(1, d)

    f_hat_val = _fhat(x_data_scaled, y_data, x_scaled, y_grid, p, q, mu, nu, bw, kernel_type)
    est = f_hat_val["est"]
    eff_n = f_hat_val["eff_n"]
    est_flag = f_hat_val["singular_flag"]

    covmat = {"cov": np.nan, "singular_flag": False}
    if cov_flag == "off":
        cov_mat = np.nan
        se = np.full(len(y_grid), np.nan)
        c_flag = False
    else:
        covmat = _cov_hat(x_data_scaled, y_data, x_scaled, y_grid, p, q, mu, nu, bw, kernel_type, cov_flag)
        cov_mat = covmat["cov"]
        c_flag = covmat["singular_flag"]
        se = np.sqrt(np.abs(np.diag(cov_mat))) * sd_y * np.mean(sd_x)

    if rbc:
        if p_RBC == p and q_RBC == q:
            est_rbc = est
            se_rbc = se
            cov_mat_rbc = cov_mat
            rbc_flag = est_flag
            c_rbc_flag = c_flag
        else:
            f_hat_rbc = _fhat(x_data_scaled, y_data, x_scaled, y_grid, p_RBC, q_RBC, mu, nu, bw, kernel_type)
            est_rbc = f_hat_rbc["est"]
            rbc_flag = f_hat_rbc["singular_flag"]
            if cov_flag == "off":
                cov_mat_rbc = np.nan
                se_rbc = np.full(len(y_grid), np.nan)
                c_rbc_flag = False
            else:
                covmat_rbc = _cov_hat(
                    x_data_scaled, y_data, x_scaled, y_grid, p_RBC, q_RBC, mu, nu, bw, kernel_type, cov_flag
                )
                cov_mat_rbc = covmat_rbc["cov"]
                c_rbc_flag = covmat_rbc["singular_flag"]
                se_rbc = np.sqrt(np.abs(np.diag(cov_mat_rbc))) * sd_y * np.mean(sd_x)
        singular_flag = bool(est_flag or c_flag or rbc_flag or c_rbc_flag)
    else:
        est_rbc = est
        se_rbc = se
        cov_mat_rbc = cov_mat
        singular_flag = bool(est_flag or c_flag)

    estimate = pd.DataFrame(
        {
            "y_grid": y_grid,
            "bw": bw,
            "est": est,
            "est_RBC": est_rbc,
            "se": se,
            "se_RBC": se_rbc,
        }
    )
    return {
        "est": estimate,
        "CovMat": {"CovMat": cov_mat, "CovMat_RBC": cov_mat_rbc},
        "x": x_scaled * sd_x + mx,
        "eff_n": eff_n,
        "singular_flag": singular_flag,
    }


def _fhat(
    x_data: np.ndarray,
    y_data: np.ndarray,
    x: np.ndarray,
    y_grid: np.ndarray,
    p: int,
    q: int,
    mu: int,
    nu: int,
    h: np.ndarray,
    kernel_type: str,
) -> dict[str, Any]:
    n = len(y_data)
    d = x_data.shape[1]
    ng = len(y_grid)
    x = x.reshape(1, d)
    singular_flag = False

    def solve_checked(mat: np.ndarray, size: int) -> np.ndarray:
        nonlocal singular_flag
        try:
            return np.linalg.inv(mat)
        except np.linalg.LinAlgError:
            singular_flag = True
            return np.zeros((size, size))

    e_nu = _basis_vec(x.reshape(-1), q, nu)
    e_mu = _basis_vec(np.array([0.0]), p, mu)
    f_hat = np.zeros(ng)
    nh_vec = np.zeros(ng)

    if np.unique(h).size == 1:
        hval = float(h[0])
        idx = np.flatnonzero(np.sum(np.abs(x_data - x) <= hval, axis=1) == d)
        x_idx = x_data[idx, :].reshape(-1, d)
        y_idx = y_data[idx]
        sort_idx = np.argsort(y_idx, kind="quicksort")
        y_sorted = y_idx[sort_idx]
        x_sorted = x_idx[sort_idx, :].reshape(-1, d)
        x_scaled = (x_sorted - x) / (hval**d)
        if x_scaled.size == 0:
            bx = np.zeros(0)
        else:
            sx_mat = solve_checked(_S_x(x_scaled, q, kernel_type) / (n * hval**d), len(e_nu))
            bx = _b_x(x_scaled, sx_mat, e_nu, q, kernel_type)

        for j, y in enumerate(y_grid):
            y_scaled = (y_sorted - y) / hval
            y_elems = np.flatnonzero(np.abs(y_scaled) <= 1)
            if len(y_elems) <= 5:
                f_hat[j] = 0.0
                nh_vec[j] = len(y_elems)
                continue
            if mu == 0:
                sx_mat = solve_checked(_S_x(x_scaled[y_elems, :].reshape(-1, d), q, kernel_type) / (n * hval**d), len(e_nu))
                bx = _b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
                sy_mat = solve_checked(_S_x(y_scaled.reshape(-1, 1), p, kernel_type) / (n * hval), len(e_mu))
                ax = _b_x(y_scaled.reshape(-1, 1), sy_mat, e_mu, p, kernel_type)
                f_hat[j] = ax[y_elems] @ np.cumsum(bx[y_elems])
            else:
                sy_mat = solve_checked(_S_x(y_scaled[y_elems].reshape(-1, 1), p, kernel_type) / (n * hval), len(e_mu))
                ax = _b_x(y_scaled[y_elems].reshape(-1, 1), sy_mat, e_mu, p, kernel_type)
                f_hat[j] = ax @ np.cumsum(bx[y_elems])
            nh_vec[j] = len(y_elems)
        f_hat = f_hat / (n**2 * hval ** (d + mu + nu + 1))
    else:
        for j, y in enumerate(y_grid):
            hval = float(h[j])
            idx = np.flatnonzero(np.sum(np.abs(x_data - x) <= hval, axis=1) == d)
            x_loc = x_data[idx, :].reshape(-1, d)
            y_loc = y_data[idx]
            sort_idx = np.argsort(y_loc, kind="quicksort")
            y_loc = y_loc[sort_idx]
            x_loc = x_loc[sort_idx, :].reshape(-1, d)
            x_scaled = (x_loc - x) / (hval**d)
            sx_mat = solve_checked(_S_x(x_scaled, q, kernel_type) / (n * hval**d), len(e_nu))
            bx = _b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
            y_scaled = (y_loc - y) / hval
            y_elems = np.flatnonzero(np.abs(y_scaled) <= 1)
            if len(y_elems) <= 5:
                f_hat[j] = 0.0
                nh_vec[j] = len(y_elems)
                continue
            if mu == 0:
                sx_mat = solve_checked(_S_x(x_scaled[y_elems, :].reshape(-1, d), q, kernel_type) / (n * hval**d), len(e_nu))
                bx = _b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
                sy_mat = solve_checked(_S_x(y_scaled.reshape(-1, 1), p, kernel_type) / (n * hval), len(e_mu))
                ax = _b_x(y_scaled.reshape(-1, 1), sy_mat, e_mu, p, kernel_type)
                f_hat[j] = ax[y_elems] @ np.cumsum(bx[y_elems])
            else:
                sy_mat = solve_checked(_S_x(y_scaled[y_elems].reshape(-1, 1), p, kernel_type) / (n * hval), len(e_mu))
                ax = _b_x(y_scaled[y_elems].reshape(-1, 1), sy_mat, e_mu, p, kernel_type)
                f_hat[j] = ax @ np.cumsum(bx[y_elems])
            f_hat[j] = f_hat[j] / (n**2 * hval ** (d + mu + nu + 1))
            nh_vec[j] = len(y_elems)

    return {"est": f_hat, "eff_n": nh_vec, "singular_flag": singular_flag}


def _cov_hat(
    x_data: np.ndarray,
    y_data: np.ndarray,
    x: np.ndarray,
    y_grid: np.ndarray,
    p: int,
    q: int,
    mu: int,
    nu: int,
    h: np.ndarray,
    kernel_type: str,
    cov_flag: str,
) -> dict[str, Any]:
    n = len(y_data)
    d = x_data.shape[1]
    ng = len(y_grid)
    x = x.reshape(1, d)
    singular_flag = False

    def solve_checked(mat: np.ndarray, size: int) -> np.ndarray:
        nonlocal singular_flag
        try:
            return np.linalg.inv(mat)
        except np.linalg.LinAlgError:
            singular_flag = True
            return np.zeros((size, size))

    e_nu = _basis_vec(x.reshape(-1), q, nu)
    e_mu = _basis_vec(np.array([1.0]), p, mu)

    def theta_at(y: float, h_val: float) -> float:
        if mu == 0:
            return float(_fhat(x_data, y_data, x, np.array([y]), 2, 1, 0, 0, np.array([h_val]), kernel_type)["est"][0])
        return float(_fhat(x_data, y_data, x, np.array([y]), p, q, mu, nu, np.array([h_val]), kernel_type)["est"][0])

    c_hat = np.zeros((ng, ng))
    if np.unique(h).size == 1:
        hval = float(h[0])
        idx = np.flatnonzero(np.sum(np.abs(x_data - x) <= hval, axis=1) == d)
        x_idx = x_data[idx, :].reshape(-1, d)
        y_idx = y_data[idx]
        sort_idx = np.argsort(y_idx, kind="quicksort")
        y_sorted = y_idx[sort_idx]
        x_sorted = x_idx[sort_idx, :].reshape(-1, d)
        x_scaled = (x_sorted - x) / (hval**d)
        if x_scaled.size == 0:
            bx = np.zeros(0)
        else:
            sx_mat = solve_checked(_S_x(x_scaled, q, kernel_type) / (n * hval**d), len(e_nu))
            bx = _b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
        theta_vals = np.array([theta_at(y, hval) for y in y_grid])
        pairs = [(i, i) for i in range(ng)] if cov_flag == "diag" else [(i, j) for i in range(ng) for j in range(i + 1)]
        if cov_flag == "diag":
            c_hat[:, :] = np.nan
        for i, j in pairs:
            y = y_grid[i]
            y_prime = y_grid[j]
            y_scaled = (y_sorted - y) / hval
            yp_scaled = (y_sorted - y_prime) / hval
            y_elems = np.flatnonzero(np.abs(y_scaled) <= 1)
            yp_elems = np.flatnonzero(np.abs(yp_scaled) <= 1)
            elems = np.intersect1d(y_elems, yp_elems, assume_unique=True)
            val = np.nan if cov_flag == "diag" else 0.0
            if len(elems) <= 5:
                val = 0.0
            elif cov_flag != "diag":
                if mu == 0:
                    sx_mat = solve_checked(_S_x((x_scaled[y_elems, :] / (n * hval**d)).reshape(-1, d), q, kernel_type), len(e_nu))
                    bx = _b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
                    sy_mat = solve_checked(_S_x(y_scaled[y_elems].reshape(-1, 1), p, kernel_type) / (n * hval), len(e_mu))
                    syp_mat = solve_checked(_S_x(yp_scaled[yp_elems].reshape(-1, 1), p, kernel_type) / (n * hval), len(e_mu))
                    a_y = _b_x(y_scaled[elems].reshape(-1, 1), sy_mat, e_mu, p, kernel_type)
                    a_yp = _b_x(yp_scaled[elems].reshape(-1, 1), syp_mat, e_mu, p, kernel_type)
                    aj = np.cumsum(a_y)
                    ak = np.cumsum(a_yp)
                else:
                    sy_mat = solve_checked(_S_x(y_scaled[y_elems].reshape(-1, 1), p, kernel_type) / (n * hval), len(e_mu))
                    syp_mat = solve_checked(_S_x(yp_scaled[yp_elems].reshape(-1, 1), p, kernel_type) / (n * hval), len(e_mu))
                    a_y = _b_x(y_scaled[y_elems].reshape(-1, 1), sy_mat, e_mu, p, kernel_type)
                    a_yp = _b_x(yp_scaled[elems].reshape(-1, 1), syp_mat, e_mu, p, kernel_type)
                    aj = np.cumsum(a_y)
                    ak = np.cumsum(a_yp)
                for k, elem in enumerate(elems):
                    kk = k + 1
                    t_1 = aj[k] * ak[k]
                    t_2 = (n - kk) * a_yp[k] * aj[k]
                    t_3 = (n - kk) * a_y[k] * ak[k]
                    t_4 = (n - kk) ** 2 * a_y[k] * a_yp[k]
                    val += bx[elem] ** 2 * (t_1 + t_2 + t_3 + t_4)
            val = val / (n * (n - 1) ** 2) - theta_vals[i] * theta_vals[j] / n**2
            c_hat[i, j] = val
            c_hat[j, i] = val
        if mu == 0:
            c_hat = c_hat * (1 / (n * hval ** (d + mu + nu))) ** 2
        else:
            c_hat = c_hat * (1 / (n * hval ** (d + mu + nu + 1))) ** 2
            np.fill_diagonal(c_hat, np.diag(c_hat / 2.0))
    else:
        theta_vals = np.array([theta_at(y_grid[i], h[i]) for i in range(ng)])
        for i in range(ng):
            for j in range(i + 1):
                hx = float(np.mean([h[i], h[j]]))
                idx = np.flatnonzero(np.sum(np.abs(x_data - x) <= hx, axis=1) == d)
                x_loc = x_data[idx, :].reshape(-1, d)
                y_loc = y_data[idx]
                sort_idx = np.argsort(y_loc, kind="quicksort")
                y_loc = y_loc[sort_idx]
                x_loc = x_loc[sort_idx, :].reshape(-1, d)
                x_scaled = (x_loc - x) / (hx**d)
                sx_mat = solve_checked(_S_x(x_scaled, q, kernel_type) / (n * hx**d), len(e_nu))
                bx = _b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
                y_scaled = (y_loc - y_grid[i]) / h[i]
                yp_scaled = (y_loc - y_grid[j]) / h[j]
                y_elems = np.flatnonzero(np.abs(y_scaled) <= 1)
                yp_elems = np.flatnonzero(np.abs(yp_scaled) <= 1)
                elems = np.intersect1d(y_elems, yp_elems, assume_unique=True)
                val = 0.0
                if len(elems) > 5:
                    if mu == 0:
                        sx_mat = solve_checked(_S_x((x_scaled[y_elems, :] / (n * hx**d)).reshape(-1, d), q, kernel_type), len(e_nu))
                        bx = _b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
                    sy_mat = solve_checked(_S_x(y_scaled[y_elems].reshape(-1, 1), p, kernel_type) / (n * h[i]), len(e_mu))
                    syp_mat = solve_checked(_S_x(yp_scaled[yp_elems].reshape(-1, 1), p, kernel_type) / (n * h[j]), len(e_mu))
                    a_y = _b_x(y_scaled[elems].reshape(-1, 1), sy_mat, e_mu, p, kernel_type)
                    a_yp = _b_x(yp_scaled[elems].reshape(-1, 1), syp_mat, e_mu, p, kernel_type)
                    aj = np.cumsum(a_y)
                    ak = np.cumsum(a_yp)
                    for k, elem in enumerate(elems):
                        kk = k + 1
                        t_1 = aj[k] * ak[k]
                        t_2 = (n - kk) * a_yp[k] * aj[k]
                        t_3 = (n - kk) * a_y[k] * ak[k]
                        t_4 = (n - kk) ** 2 * a_y[k] * a_yp[k]
                        val += bx[elem] ** 2 * (t_1 + t_2 + t_3 + t_4)
                val = val / (n * (n - 1) ** 2) - theta_vals[i] * theta_vals[j] / n**2
                c_hat[i, j] = val
                c_hat[j, i] = val
        if mu == 0:
            scale = 1 / (n * h ** (d + mu + nu))
        else:
            scale = 1 / (n * h ** (d + mu + nu + 1))
        c_hat = c_hat * scale[:, None] * scale[None, :]
        if mu != 0:
            np.fill_diagonal(c_hat, np.diag(c_hat / 2.0))

    return {"cov": c_hat, "singular_flag": singular_flag}


def _bw_rot(
    y_data: np.ndarray,
    x_data: np.ndarray,
    y_grid: np.ndarray,
    x: np.ndarray,
    p: int,
    q: int,
    mu: int,
    nu: int,
    kernel_type: str,
    regularize: bool,
) -> np.ndarray:
    return _bw_rot_common(y_data, x_data, y_grid, x, p, q, mu, nu, kernel_type, regularize, integrated=False)


def _bw_irot(
    y_data: np.ndarray,
    x_data: np.ndarray,
    y_grid: np.ndarray,
    x: np.ndarray,
    p: int,
    q: int,
    mu: int,
    nu: int,
    kernel_type: str,
    regularize: bool,
) -> np.ndarray:
    return _bw_rot_common(y_data, x_data, y_grid, x, p, q, mu, nu, kernel_type, regularize, integrated=True)


def _bw_rot_common(
    y_data: np.ndarray,
    x_data: np.ndarray,
    y_grid: np.ndarray,
    x: np.ndarray,
    p: int,
    q: int,
    mu: int,
    nu: int,
    kernel_type: str,
    regularize: bool,
    integrated: bool,
) -> np.ndarray:
    sd_y = _sd(y_data)
    sd_x = np.apply_along_axis(_sd, 0, x_data)
    mx = np.mean(x_data, axis=0)
    my = float(np.mean(y_data))
    d = x_data.shape[1]
    n = len(y_data)
    ng = len(y_grid)
    data = np.column_stack([y_data, x_data])
    x = _as_vector(x)
    if d == 1:
        bx_band = 0.5 if integrated else 1.06 * n ** (-1 / 5) * sd_x[0]
        lower_x = float(np.min(x_data[:, 0]) - x[0])
        upper_x = float(np.max(x_data[:, 0]) - x[0])
        bias_dgp = np.full((ng, 3), np.nan)
        for j, y in enumerate(y_grid):
            lower_y = float(np.min(y_data) - y)
            upper_y = float(np.max(y_data) - y)
            bias_dgp[j, 0] = _normal_dgps(y, mu, my, sd_y) * _normal_dgps(x[0], 2, mx[0], sd_x[0])
            bias_dgp[j, 1] = _normal_dgps(y, p + 1, my, sd_y) * _normal_dgps(x[0], 0, mx[0], sd_x[0])
            Sx = np.linalg.inv(_S_exact(lower_x, upper_x, np.array([x[0]]), q, kernel_type))
            Sy = np.linalg.inv(_S_exact(lower_y, upper_y, np.array([y]), p, kernel_type))
            cx = _c_exact(lower_x, upper_x, np.array([x[0]]), q + 1, q, kernel_type)
            cy = _c_exact(lower_y, upper_y, np.array([y]), p + 1, p, kernel_type)
            bias_dgp[j, 0] *= (cx.T @ Sx).reshape(-1)[nu]
            bias_dgp[j, 1] *= (cy.T @ Sy).reshape(-1)[mu]
            bias_dgp[j, 2] = (bias_dgp[j, 0] + bias_dgp[j, 1]) ** 2

        v_dgp = np.full(ng, np.nan)
        for j, y in enumerate(y_grid):
            lower_y = float(np.min(y_data) - y)
            upper_y = float(np.max(y_data) - y)
            Sx = np.linalg.inv(_S_exact(lower_x, upper_x, np.array([x[0]]), q, kernel_type))
            if mu > 0:
                Sy = np.linalg.inv(_S_exact(lower_y, upper_y, np.array([y]), p, kernel_type))
                Ty = _T_y_exact(lower_y, upper_y, np.array([y]), p)
                if mu == 1:
                    Tx = _T_x(x_data, np.array([x[0]]), q, np.array([bx_band]), kernel_type) / bx_band
                else:
                    Tx = _T_y_exact(lower_x, upper_x, np.array([x[0]]), q)
                v_dgp[j] = norm.pdf(y) * norm.cdf(x[0]) * (Sy @ Ty @ Sy)[mu, mu] * (Sx @ Tx @ Sx)[nu, nu]
            else:
                Tx = _T_x(x_data, np.array([x[0]]), q, np.array([bx_band]), kernel_type) / bx_band
                cdf_hat = norm.cdf(y) * norm.cdf(x[0])
                v_dgp[j] = cdf_hat * (1 - cdf_hat) * (Sx @ Tx @ Sx)[nu, nu]

        if mu == 0:
            alpha = d + 2 * min(p, q) + 2 * min(mu, nu) + 2 * nu + 2
        else:
            alpha = d + 2 * min(p, q) + 2 * max(mu, nu) + 1
        if integrated:
            h = (abs(np.mean(v_dgp) / (2 * np.mean(bias_dgp[:, 2])))) ** (1 / alpha) * n ** (-1 / alpha)
            h = sd_y * sd_x[0] * h
            if regularize:
                h = max(h, np.sort(np.abs(x_data[:, 0] - x[0]))[min(n, 20 + q + 1) - 1])
                for y in y_grid:
                    h = max(h, np.sort(np.abs(y_data - y))[min(n, 20 + p + 1) - 1])
            return np.array([h])
        h = (np.abs(v_dgp / bias_dgp[:, 2])) ** (1 / alpha) * n ** (-1 / alpha)
        if regularize:
            for j, y in enumerate(y_grid):
                h[j] = max(h[j], np.sort(np.abs(x_data[:, 0] - x[0]))[min(n, 20 + q + 4) - 1])
                h[j] = max(h[j], np.sort(np.abs(y_data - y))[min(n, 20 + p + 4) - 1])
        return h

    bx_band = (4 / (d + 1)) ** (1 / (d + 4)) * n ** (-1 / (d + 4)) * sd_x
    sigma_hat = np.cov(x_data, rowvar=False, ddof=1)
    sigma_inv = np.linalg.inv(sigma_hat)
    x_minus_mx = (x - mx).reshape(d, 1)
    mvt_deriv = _dmvnorm(x, np.zeros(d), sigma_hat) * (sigma_inv @ x_minus_mx @ x_minus_mx.T @ sigma_inv - sigma_inv)
    corr_x = np.corrcoef(x_data, rowvar=False)
    cov_data = np.cov(data, rowvar=False, ddof=1)
    bias_dgp = np.zeros((ng, 3))
    for j, y in enumerate(y_grid):
        bias_dgp[j, 1] = _normal_dgps(y, 3, my, sd_y) * _pmvnorm(x, np.zeros(d), corr_x)
        Sx = np.linalg.inv(_S_exact(eval_pt=x, p=q, kernel_type=kernel_type))
        Sy = np.linalg.inv(_S_exact(eval_pt=np.array([y]), p=p, kernel_type=kernel_type))
        cy = _c_exact(eval_pt=np.array([y]), m=p + 1, p=p, kernel_type=kernel_type)
        b1 = np.exp(-0.5 * (y - my) ** 2 / sd_y) / (np.sqrt(2 * np.pi) * sd_y) * (y - my) / sd_y * mvt_deriv
        diag_elems = np.diag(b1).copy()
        upper_vals = b1.copy()
        upper_vals[np.tril_indices(d)] = np.nan
        pxy = upper_vals.reshape(-1, order="F")[~np.isnan(upper_vals.reshape(-1, order="F"))]
        mv = _mvec(q + 1, d)
        for i in range(len(mv) - d):
            cx = _c_x(x_data, x, np.asarray(mv[i]), q, bx_band, kernel_type)
            bias_dgp[j, 0] += pxy[i] * (cx.T @ Sx).reshape(-1)[nu]
        for i in range(d):
            cx = _c_x(x_data, x, np.asarray(mv[len(mv) - d + i]), q, bx_band, kernel_type)
            bias_dgp[j, 0] += diag_elems[i] * (cx.T @ Sx).reshape(-1)[nu]
        bias_dgp[j, 1] *= (cy.T @ Sy).reshape(-1)[mu]
        if integrated:
            z = np.r_[y, x]
            bias_dgp[j, 2] = _dmvnorm(z, np.zeros(d + 1), cov_data) * (bias_dgp[j, 0] + bias_dgp[j, 1]) ** 2
        else:
            bias_dgp[j, 2] = (bias_dgp[j, 0] + bias_dgp[j, 1]) ** 2

    v_dgp = np.full(ng, np.nan)
    Sx = np.linalg.inv(_S_exact(eval_pt=x, p=q, kernel_type=kernel_type))
    for j, y in enumerate(y_grid):
        Sy = np.linalg.inv(_S_exact(eval_pt=np.array([y]), p=p, kernel_type=kernel_type))
        Ty = _T_y_exact(eval_pt=np.array([y]), p=p)
        Tx = _T_x(x_data, x, q, bx_band, kernel_type)
        z = np.r_[y, x]
        if mu > 0:
            v_dgp[j] = _dmvnorm(z, np.zeros(d + 1), cov_data) ** 2
            v_dgp[j] *= (Sy @ Ty @ Sy)[mu, mu] * (Sx @ Tx @ Sx)[nu, nu]
        else:
            cdf_hat = _pmvnorm(z, np.zeros(d + 1), cov_data)
            v_dgp[j] = cdf_hat * (1 - cdf_hat) * _dmvnorm(z, np.zeros(d + 1), cov_data)
            v_dgp[j] *= (Sx @ Tx @ Sx)[nu, nu]

    if integrated:
        h = (np.sum(v_dgp) / np.sum(bias_dgp[:, 2])) ** (1 / 6) * n ** (-1 / 6)
    else:
        h = (v_dgp / bias_dgp[:, 2]) ** (1 / 6) * n ** (-1 / 6)
    h = _sd(y_data) * _sd(x_data.reshape(-1, order="F")) * h
    return np.asarray(h)


def _S_x(x_data: np.ndarray, q: int, kernel_type: str) -> np.ndarray:
    x_data = _as_matrix(x_data)
    if x_data.shape[1] == 1:
        x = x_data[:, 0]
        r = _univar_basis_matrix(x, q)
        k = _kernel_eval_vec(x, kernel_type)
        return r.T @ (r * k[:, None])
    r = np.column_stack([_poly_base(row, q) for row in x_data])
    k = np.array([_kernel_eval(row, kernel_type) for row in x_data])
    return r @ (r.T * k[:, None])


def _T_x(x_data: np.ndarray, eval_pt: np.ndarray, q: int, h: np.ndarray, kernel_type: str) -> np.ndarray:
    x_data = _as_matrix(x_data)
    eval_pt = _as_vector(eval_pt)
    h = _as_vector(h)
    d = x_data.shape[1]
    n = x_data.shape[0]
    e_base = _basis_vec(eval_pt, q, 0)
    if d == 1:
        x_pol = (x_data[:, 0] - eval_pt[0]) / h[0]
        r = _univar_basis_matrix(x_pol, q)
        k = _kernel_eval_vec(x_pol, kernel_type)
        rk = r * k[:, None]
        return rk.T @ rk / n
    x_pol = _r_recycle_binary(x_data - eval_pt[0], h, op="/")
    r = np.column_stack([_poly_base(row, q) for row in x_pol])
    k = np.array([_kernel_eval(row, kernel_type) for row in x_pol])
    return ((r.T * k[:, None]).T @ (r.T * k[:, None])) / n


def _T_y_exact(
    lower: float = -1,
    upper: float = 1,
    eval_pt: np.ndarray | None = None,
    p: int = 1,
    kernel_type: str = "uniform",
) -> np.ndarray:
    v = np.zeros((p + 1, p + 1))
    a = max(lower, -1)
    b = min(upper, 1)
    if a < b:
        for i in range(p + 1):
            for j in range(p + 1):
                if kernel_type == "uniform":
                    num1 = -(b ** (i + j + 3) - a ** (i + j + 3))
                    denom1 = (i + 1) * (i + 2) * (i + j + 3)
                    num2 = -(a ** (i + 2) * (b ** (j + 1) - a ** (j + 1)))
                    denom2 = (i + 2) * (j + 1)
                    num3 = b ** (i + 1) * (b ** (j + 2) - a ** (j + 2))
                    denom3 = (i + 2) * (j + 2)
                    v[i, j] = (num1 / denom1 + num2 / denom2 + num3 / denom3) / 4
    fact = np.outer([factorial(i) for i in range(p + 1)], [factorial(i) for i in range(p + 1)])
    return v / fact


def _S_exact(
    lower: float = -1,
    upper: float = 1,
    eval_pt: np.ndarray | None = None,
    p: int = 1,
    kernel_type: str = "epanechnikov",
) -> np.ndarray:
    if eval_pt is None:
        eval_pt = np.array([0.0])
    eval_pt = _as_vector(eval_pt)
    if len(eval_pt) == 1:
        lower_lim = max(lower, -1)
        upper_lim = min(upper, 1)
        poly_list = np.zeros(2 * p + 1)
        if lower_lim < upper_lim:
            for i in range(2 * p + 1):
                poly_list[i] = _int_val(i, lower_lim, upper_lim, kernel_type)
        poly_mat = np.zeros((p + 1, p + 1))
        for j in range(p + 1):
            poly_mat[:, j] = poly_list[j : j + p + 1]
        fact = np.outer([factorial(i) for i in range(p + 1)], [factorial(i) for i in range(p + 1)])
        return norm.pdf(eval_pt[0]) * poly_mat / fact
    poly_list = [0]
    for i in range(1, p + 1):
        poly_list.extend([i] * int(np.sum(_basis_vec(eval_pt, p, i))))
    poly_mat = np.zeros((len(poly_list), len(poly_list)))
    poly_mat[:, 0] = poly_list
    poly_mat[0, :] = poly_list
    for i in range(1, len(poly_list)):
        for j in range(1, len(poly_list)):
            val = _int_val(poly_mat[0, j] + poly_mat[i, 0], -1, 1, kernel_type)
            denom = factorial(int(poly_mat[0, j])) * factorial(int(poly_mat[i, 0]))
            poly_mat[i, j] = val / denom if i == j else val**2 / denom
    return poly_mat


def _c_exact(
    lower: float = -1,
    upper: float = -1,
    eval_pt: np.ndarray | None = None,
    m: int = 1,
    p: int = 1,
    kernel_type: str = "epanechnikov",
) -> np.ndarray:
    if eval_pt is None:
        eval_pt = np.array([0.0])
    eval_pt = _as_vector(eval_pt)
    if len(eval_pt) == 1:
        v = np.zeros((p + 1, 1))
        lower_lim = max(lower, -1)
        upper_lim = min(upper, 1)
        if lower_lim < upper_lim:
            for i in range(p + 1):
                v[i, 0] = _int_val(i + m, lower_lim, upper_lim, kernel_type) / factorial(m)
        return v
    return np.array([[np.sum(_basis_vec(eval_pt, p, m))]])


def _c_x(x_data: np.ndarray, eval_pt: np.ndarray, m: np.ndarray, q: int, h: np.ndarray, kernel_type: str) -> np.ndarray:
    x_data = _as_matrix(x_data)
    eval_pt = _as_vector(eval_pt)
    h = _as_vector(h)
    n = x_data.shape[0]
    c = np.zeros((1, len(_basis_vec(eval_pt, q, 0))))
    v = _r_recycle_binary(x_data - _r_recycle_matrix(eval_pt, x_data.shape), h, op="/")
    x_pol = x_data - _r_recycle_matrix(eval_pt / h, x_data.shape)
    r = np.column_stack([_poly_base(row, q) for row in v])
    m_mat = _r_recycle_matrix(m, x_data.shape)
    fact_mat = _r_recycle_matrix(np.array([factorial(int(vv)) for vv in m]), x_data.shape)
    m_poly = np.sum((x_pol**m_mat) / fact_mat, axis=1)
    m_poly = _r_recycle_binary(m_poly, h, op="/") * x_data.shape[1]
    k = np.array([_kernel_eval(row, kernel_type) for row in v])
    return (r.T * k[:, None]).T @ m_poly.reshape(-1, 1) / n


def _b_x(datavec: np.ndarray, s_mat: np.ndarray, e_vec: np.ndarray, q: int, kernel_type: str) -> np.ndarray:
    datavec = _as_matrix(datavec)
    if datavec.shape[0] == 0:
        return np.zeros(0)
    if datavec.shape[1] == 1:
        x = datavec[:, 0]
        r = _univar_basis_matrix(x, q)
        k = _kernel_eval_vec(x, kernel_type)
        return ((r * k[:, None]) @ (s_mat @ e_vec)).reshape(-1)
    weights = s_mat @ e_vec
    out = np.zeros(datavec.shape[0])
    for i, row in enumerate(datavec):
        out[i] = (_poly_base(row, q) * _kernel_eval(row, kernel_type)) @ weights
    return out


def _poly_base(x: ArrayLike, p: int) -> np.ndarray:
    x = _as_vector(x)
    if len(x) == 1:
        return np.array([x[0] ** i / factorial(i) for i in range(p + 1)], dtype=float)
    if p >= 5:
        raise ValueError("Not implementable for multivariate polynomials of order greater than 4.")
    polyvec: list[float] = [1.0]
    powers_by_var = [x ** i for i in range(1, p + 1)]
    if p >= 1:
        polyvec.extend(powers_by_var[0].tolist())
    if p >= 2:
        polyvec.extend((powers_by_var[1] / factorial(2)).tolist())
        for m in range(2, p + 1):
            for powers in _print_all_sum(m):
                values = []
                counts = {power: powers.count(power) for power in set(powers)}
                repeated_parts: list[list[float]] = []
                keep_parts: list[np.ndarray] = []
                for power, count in counts.items():
                    if count > 1:
                        arr = powers_by_var[power - 1]
                        repeated_parts.append([float(np.prod(arr[list(idx)])) for idx in _combinations(range(len(arr)), count)])
                    else:
                        keep_parts.append(powers_by_var[power - 1])
                if keep_parts:
                    for vals in product(*[part.tolist() for part in keep_parts], *repeated_parts):
                        values.append(float(np.prod(vals)))
                else:
                    values = [v for part in repeated_parts for v in part]
                polyvec.extend([v / factorial(m) for v in values])
    return np.asarray(polyvec, dtype=float)


def _basis_vec(x: ArrayLike, p: int, mu: int) -> np.ndarray:
    x = _as_vector(x)
    if len(x) == 1:
        out = np.zeros(p + 1)
        out[mu] = 1.0
        return out
    # The R implementation marks all entries of the requested total degree.
    lengths = [1]
    if p >= 1:
        lengths.append(len(x))
    if p >= 2:
        lengths.append(len(x))
        for m in range(2, p + 1):
            for powers in _print_all_sum(m):
                count = 1
                counts = {power: powers.count(power) for power in set(powers)}
                repeated_sizes = []
                keep_count = 0
                for power, freq in counts.items():
                    if freq > 1:
                        repeated_sizes.append(len(list(_combinations(range(len(x)), freq))))
                    else:
                        keep_count += 1
                size = int(np.prod([len(x)] * keep_count + repeated_sizes)) if keep_count or repeated_sizes else 0
                lengths.append(size)
    out = np.zeros(sum(lengths))
    cursor = 0
    for degree, size in enumerate(lengths):
        if degree == mu:
            out[cursor : cursor + size] = 1.0
        cursor += size
    return out


def _univar_basis_matrix(x: np.ndarray, q: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    r = np.ones((len(x), q + 1))
    for j in range(1, q + 1):
        r[:, j] = r[:, j - 1] * x / j
    return r


def _kernel_eval(x: ArrayLike, kernel_type: str) -> float:
    x = np.asarray(x, dtype=float)
    if kernel_type == "uniform":
        k = np.where(np.abs(x) <= 1, 0.5, 0.0)
    elif kernel_type == "triangular":
        k = (1 - np.abs(x)) * np.where(np.abs(x) <= 1, 1.0, 0.0)
    else:
        k = 0.75 * (1 - x**2) * np.where(np.abs(x) <= 1, 1.0, 0.0)
    return float(np.prod(k))


def _kernel_eval_vec(x: np.ndarray, kernel_type: str) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if kernel_type == "uniform":
        return np.where(np.abs(x) <= 1, 0.5, 0.0)
    if kernel_type == "triangular":
        return (1 - np.abs(x)) * np.where(np.abs(x) <= 1, 1.0, 0.0)
    return 0.75 * (1 - x**2) * np.where(np.abs(x) <= 1, 1.0, 0.0)


def _int_val(l: float, a: float, b: float, kernel_type: str) -> float:
    l = int(l)
    if kernel_type == "triangular":
        if a >= 0 and b >= 0:
            num = a ** (l + 1) * (-2 + a + (-1 + a) * l) + b ** (l + 1) * (2 + l - b * (1 + l))
            return num / ((l + 1) * (l + 2))
        if a < 0 and b > 0:
            num1 = -a ** (l + 1) * (2 + a + l + a * l)
            denom1 = 2 + 3 * l + l * l
            num2 = (2 + l) * b ** (l + 1) - (l + 1) * b ** (l + 2)
            denom2 = (l + 1) * (l + 2)
            return num1 / denom1 + num2 / denom2
        if a < 0 and b < 0:
            num = -a ** (l + 1) * (2 + a + l + a * l) + b ** (l + 1) * (2 + l + b + b * l)
            return num / ((l + 1) * (l + 2))
    if kernel_type == "uniform":
        return 0.5 * (b ** (l + 1) - a ** (l + 1)) / (l + 1)
    num = (l + 1) * a ** (l + 3) - (l + 3) * a ** (l + 1) + b ** (l + 1) * (3 + l - (b**2) * (l + 1))
    return 0.75 * num / ((l + 1) * (l + 3))


def _normal_dgps(x: float, v: int, mean: float, sd: float) -> float:
    if v == 0:
        return float(norm.cdf(x, loc=mean, scale=sd))
    z = (x - mean) / sd
    pdf = norm.pdf(z) / sd
    order = v - 1
    return float(((-1) ** order) * eval_hermitenorm(order, z) * pdf / (sd**order))


def _dmvnorm(x: np.ndarray, mean: np.ndarray, sigma: np.ndarray) -> float:
    return float(multivariate_normal.pdf(np.asarray(x, dtype=float), mean=mean, cov=sigma, allow_singular=False))


def _pmvnorm(x: np.ndarray, mean: np.ndarray, sigma: np.ndarray) -> float:
    return float(multivariate_normal.cdf(np.asarray(x, dtype=float), mean=mean, cov=sigma, allow_singular=False))


def _mvec(n: int, d: int) -> list[list[int]]:
    if d == 1:
        return [[n]]
    pvec = _print_all_sum(n)
    padded = [vals + [0] * (d - len(vals)) for vals in pvec]
    out: list[tuple[int, ...]] = []
    for vals in padded:
        out.extend(permutations(vals))
    out.extend(permutations([n] + [0] * (d - 1)))
    unique: list[tuple[int, ...]] = []
    for vals in out:
        if len(vals) <= d and vals not in unique:
            unique.append(vals)
    return [list(vals) for vals in unique]


def _print_all_sum(target: int) -> list[list[int]]:
    output: list[list[int]] = []

    def rec(current_sum: int, start: int, result: list[int]) -> None:
        if target == current_sum:
            output.append(result.copy())
        for i in range(start, target):
            temp_sum = current_sum + i
            if temp_sum <= target:
                result.append(i)
                rec(temp_sum, i, result)
                result.pop()
            else:
                return

    rec(0, 1, [])
    return output


def _combinations(values: range, r: int) -> list[tuple[int, ...]]:
    return list(combinations(values, r))


def _r_recycle_matrix(values: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    flat_len = shape[0] * shape[1]
    recycled = np.resize(np.asarray(values, dtype=float), flat_len)
    return recycled.reshape(shape, order="F")


def _r_recycle_binary(left: np.ndarray, right: np.ndarray, op: str) -> np.ndarray:
    left_arr = np.asarray(left, dtype=float)
    right_arr = np.resize(np.asarray(right, dtype=float), left_arr.size).reshape(left_arr.shape, order="F")
    if op == "/":
        return left_arr / right_arr
    if op == "*":
        return left_arr * right_arr
    if op == "-":
        return left_arr - right_arr
    raise ValueError(op)


def _r_quantile(x: np.ndarray, probs: ArrayLike) -> np.ndarray:
    return np.quantile(np.asarray(x, dtype=float), probs, method="linear")


def _as_vector(x: ArrayLike) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    return arr.reshape(-1)


def _as_matrix(x: ArrayLike) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1, 1)
    elif arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    elif arr.ndim > 2:
        raise ValueError("Expected a vector or matrix.")
    return arr


def _sd(x: np.ndarray) -> float:
    return float(np.std(np.asarray(x, dtype=float), ddof=1))


def _check_order(name: str, value: int) -> None:
    if not isinstance(value, (int, np.integer)) or value < 0 or value > 20:
        raise ValueError(f"Polynomial or derivative order {name} incorrectly specified.")


def _check_kernel(kernel_type: str) -> str:
    kernel_type = (kernel_type or "epanechnikov").lower()
    if kernel_type not in {"triangular", "uniform", "epanechnikov"}:
        raise ValueError("Kernel function incorrectly specified.")
    return kernel_type


def load_r_fixture_rows(path: str | Path) -> pd.DataFrame:
    out = pd.read_csv(path)
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    return out
