# -------------------------------------------------------------------
# LPCDE Python replication illustration
# -------------------------------------------------------------------

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from lpcde import confint, lpbwcde, lpcde, plot, summary


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "Python" / "lpcde" / "tests" / "fixtures" / "r-baseline-1d.csv"


def main() -> None:
    data = pd.read_csv(FIXTURE)
    x_data = data[["x"]].to_numpy(dtype=float)
    y_data = data["y"].to_numpy(dtype=float)
    y_grid = np.arange(-2, 2.1, 1, dtype=float)

    model1 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1)
    print(summary(model1))
    print(confint(model1))

    model_reg = lpcde(
        x_data=x_data,
        y_data=y_data,
        y_grid=y_grid,
        x=0,
        bw=1,
        nonneg=True,
        normalize=True,
    )
    print(summary(model_reg))

    model_deriv = lpcde(
        x_data=x_data,
        y_data=y_data,
        y_grid=y_grid,
        x=0,
        p=3,
        mu=2,
        bw=1,
    )
    print(summary(model_deriv))

    yq = np.quantile(y_data, np.arange(0.2, 0.81, 0.2), method="linear")
    bw_mse = lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid=yq)
    print(summary(bw_mse))

    bw_imse = lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid=yq, bw_type="imse-rot")
    print(summary(bw_imse))

    plot(model1)


if __name__ == "__main__":
    main()
