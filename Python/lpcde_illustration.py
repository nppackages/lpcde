#--------------------------------------------------------------------------------
# LPCDE Package
# Illustration Code
#--------------------------------------------------------------------------------

import numpy as np
import matplotlib

matplotlib.use("Agg")

from lpcde import confint, lpbwcde, lpcde, plot, summary


#----------------------------------------
# Generate data
#----------------------------------------
np.random.seed(42)
n = 500
x_data = np.random.normal(size=n)
y_data = np.random.normal(loc=x_data, scale=1)
y_grid = np.linspace(-1, 1, 5)


#-------------------------------------------------------------------
# lpcde(): estimation with bandwidth 0.5 on provided grid points
#-------------------------------------------------------------------
model1 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=0.5)
print(summary(model1))
print(confint(model1))


#----------------------------------------
# lpbwcde(): bandwidth selection
#----------------------------------------
model1bw = lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid=y_grid)
print(summary(model1bw))


#----------------------------------------
# plot(): base matplotlib output
#----------------------------------------
plot(model1)
