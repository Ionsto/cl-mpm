import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import os
import numpy as np
import re

from scipy import integrate

top_dir = "../../../"
regex = re.compile(r'^output.*')
folders = list(filter(regex.search,os.listdir(top_dir)))
plt.figure(1)
def get_load(filename):
    mpm = pd.read_csv(top_dir+filename)
    mpm["time"] = mpm["time"]
    mpm["disp"] = mpm["disp"]
    return mpm


folders.sort()

for i in folders:
    print("loading folder: ",i)
    mpm = get_load("./{}/disp.csv".format(i))
    if len(mpm["time"]) > 0:
        plt.plot(mpm["time"].values,mpm["disp"].values,label=i)
plt.xlabel("Displacement (mm)")
plt.ylabel("Load (N)")
plt.legend()
plt.show()
