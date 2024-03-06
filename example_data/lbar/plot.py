import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import re

from scipy import integrate
def calculate_gf(disp,load):
    return integrate.trapz(load,disp)#/(0.102*0.6*13e-3)

data = pd.read_csv("load-disp.csv")


print("GF experimental:",calculate_gf(1e0*data["disp"],100e-3*data["load"]))


regex = re.compile(r'output-.*')
folders = list(filter(regex.search,os.listdir("./")))

plt.figure()
plt.plot(data["disp"],data["load"],label="Data")

data = pd.read_csv("mvm.csv")
plt.plot(data["x"],1e3*data["y"],label="Modified VM")

lower = pd.read_csv("lower.csv")
upper = pd.read_csv("upper.csv")
x_min = min(lower["x"].min(),upper["x"].min())
x_max = max(lower["x"].max(),upper["x"].max())
x_samples = np.linspace(x_min,x_max,100)

lower_y = np.interp(x_samples,lower["x"],1e3*lower["y"])
upper_y = np.interp(x_samples,upper["x"],1e3*upper["y"])

plt.fill_between(x_samples,lower_y,upper_y)

for i in folders:
    print("loading folder: ",i)
    mpm = pd.read_csv("./{}/disp.csv".format(i))
    print("GF mpm:",calculate_gf(1e3*mpm["disp"],100e-3*mpm["load"]))
    plt.plot(1e3*mpm["disp"],100e-3*mpm["load"],label=i)
    plt.xlabel("Displacement (mm)")
    plt.ylabel("Load (N)")
    plt.legend()
plt.show()
