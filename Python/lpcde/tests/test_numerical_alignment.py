from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lpcde import lpbwcde, lpcde


FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _mk_matrix(case, prefix, value):
    if isinstance(value, pd.DataFrame):
        arr = value.to_numpy(dtype=float)
        cols = list(value.columns)
    else:
        arr = np.asarray(value, dtype=float)
        cols = list(range(1, arr.shape[1] + 1))
    rows = []
    for col_idx, col in enumerate(cols):
        for row_idx in range(arr.shape[0]):
            rows.append(
                {
                    "case": case,
                    "metric": f"{prefix}[{row_idx + 1},{col}]",
                    "value": arr[row_idx, col_idx],
                }
            )
    return rows


def _mk_value(case, prefix, value):
    if np.isscalar(value):
        return [{"case": case, "metric": prefix, "value": value}]
    return _mk_matrix(case, prefix, value)


def _mk_bw(case, value):
    return _mk_matrix(case, "BW", value)


def _make_numerical_baseline():
    one = pd.read_csv(FIXTURES / "r-baseline-1d.csv")
    x_data = one[["x"]].to_numpy()
    y_data = one["y"].to_numpy()
    y_grid = np.arange(-2, 3, 1, dtype=float)

    model1 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1)
    model_reg = lpcde(
        x_data=x_data,
        y_data=y_data,
        y_grid=y_grid,
        x=0,
        bw=1,
        nonneg=True,
        normalize=True,
    )
    model_deriv = lpcde(
        x_data=x_data,
        y_data=y_data,
        y_grid=y_grid,
        x=0,
        p=3,
        mu=2,
        bw=1,
    )
    model_tri_diag = lpcde(
        x_data=x_data,
        y_data=y_data,
        y_grid=y_grid,
        x=0,
        bw=1.1,
        kernel_type="triangular",
        cov_flag="diag",
    )
    model_unif_off = lpcde(
        x_data=x_data,
        y_data=y_data,
        y_grid=y_grid,
        x=0,
        bw=1.2,
        kernel_type="uniform",
        cov_flag="off",
    )

    yq = np.quantile(y_data, np.arange(0.2, 0.81, 0.2), method="linear")
    bw_mse = lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid=yq)
    bw_imse = lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid=yq, bw_type="imse-rot")

    two = pd.read_csv(FIXTURES / "r-baseline-2d.csv")
    x2_data = two[["x1", "x2"]].to_numpy()
    y2_data = two["y"].to_numpy()
    y2_grid = np.quantile(y2_data, [0.25, 0.5, 0.75], method="linear")
    model_2d = lpcde(
        x_data=x2_data,
        y_data=y2_data,
        y_grid=y2_grid,
        x=np.array([[0, 0]]),
        bw=1,
        cov_flag="diag",
    )
    bw_2d = lpbwcde(
        y_data=y2_data,
        x_data=x2_data,
        x=np.array([[0, 0]]),
        y_grid=y2_grid,
        bw_type="imse-rot",
    )

    rows = []
    rows += _mk_matrix("fixed_bw", "Estimate", model1.estimate)
    rows += _mk_matrix("fixed_bw", "CovMat", model1.cov_mat["CovMat"])
    rows += _mk_matrix("fixed_bw", "CovMat_RBC", model1.cov_mat["CovMat_RBC"])
    rows += _mk_matrix("regularized_bw", "Estimate", model_reg.estimate)
    rows += _mk_matrix("regularized_bw", "CovMat", model_reg.cov_mat["CovMat"])
    rows += _mk_matrix("regularized_bw", "CovMat_RBC", model_reg.cov_mat["CovMat_RBC"])
    rows += _mk_matrix("derivative_mu2", "Estimate", model_deriv.estimate)
    rows += _mk_matrix("derivative_mu2", "CovMat", model_deriv.cov_mat["CovMat"])
    rows += _mk_matrix("derivative_mu2", "CovMat_RBC", model_deriv.cov_mat["CovMat_RBC"])
    rows += _mk_matrix("triangular_diag", "Estimate", model_tri_diag.estimate)
    rows += _mk_matrix("triangular_diag", "CovMat", model_tri_diag.cov_mat["CovMat"])
    rows += _mk_matrix("triangular_diag", "CovMat_RBC", model_tri_diag.cov_mat["CovMat_RBC"])
    rows += _mk_matrix("uniform_off", "Estimate", model_unif_off.estimate)
    rows += _mk_value("uniform_off", "CovMat", model_unif_off.cov_mat["CovMat"])
    rows += _mk_value("uniform_off", "CovMat_RBC", model_unif_off.cov_mat["CovMat_RBC"])
    rows += _mk_matrix("two_dimensional_diag", "Estimate", model_2d.estimate)
    rows += _mk_matrix("two_dimensional_diag", "CovMat", model_2d.cov_mat["CovMat"])
    rows += _mk_matrix("two_dimensional_diag", "CovMat_RBC", model_2d.cov_mat["CovMat_RBC"])
    rows += _mk_bw("bandwidth_mse", bw_mse.bw)
    rows += _mk_bw("bandwidth_imse", bw_imse.bw)
    rows += _mk_bw("bandwidth_2d_imse", bw_2d.bw)
    return pd.DataFrame(rows)


def test_python_matches_r_numerical_baseline():
    fixture = pd.read_csv(FIXTURES / "r-baseline-output.csv")
    fixture["value"] = pd.to_numeric(fixture["value"], errors="coerce")
    current = _make_numerical_baseline()

    fixture_key = fixture["case"] + "::" + fixture["metric"]
    current_key = current["case"] + "::" + current["metric"]
    assert set(current_key) == set(fixture_key)

    current = current.set_index(current_key).loc[fixture_key].reset_index(drop=True)
    expected = fixture["value"].to_numpy()
    actual = current["value"].to_numpy()
    ok = (np.isnan(expected) & np.isnan(actual)) | np.isclose(actual, expected, rtol=0, atol=1e-8)
    if not np.all(ok):
        bad = np.flatnonzero(~ok)[:10]
        detail = [
            f"{fixture_key.iloc[i]}: actual={actual[i]:.17g}, expected={expected[i]:.17g}"
            for i in bad
        ]
        pytest.fail("Python/R numerical mismatch:\n" + "\n".join(detail))
