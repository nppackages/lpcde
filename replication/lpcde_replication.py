# -------------------------------------------------------------------
# LPCDE Python replication script for the software article examples
# -------------------------------------------------------------------

from pathlib import Path

import numpy as np
import pandas as pd

from lpcde import lpbwcde, lpcde, summary


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "Python" / "lpcde" / "tests" / "fixtures" / "r-baseline-1d.csv"


def main() -> None:
    data = pd.read_csv(FIXTURE)
    x_data = data[["x"]].to_numpy(dtype=float)
    y_data = data["y"].to_numpy(dtype=float)

    y_grid = np.arange(-2, 2.1, 0.5, dtype=float)
    model1 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1, rbc=True)
    print(summary(model1))

    model_reg = lpcde(
        x_data=x_data,
        y_data=y_data,
        y_grid=y_grid,
        x=0,
        bw=1,
        rbc=True,
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

    bw = lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid=y_grid)
    print(summary(bw))


if __name__ == "__main__":
    main()
