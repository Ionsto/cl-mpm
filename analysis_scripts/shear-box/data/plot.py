import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

load_zeroing = True
# load_zeroing = False
# load_combined = True
load_combined = False
load_clipping = False

def get_load(filename):
    mpm = pd.read_csv(filename)
    if load_clipping:
        mpm = mpm[mpm["disp"] >= 0.01e-3]
    if len(mpm["load"]) > 0:
        mpm["load-diff"] = mpm["l-left"] + mpm["l-right"]
        if load_combined:
            mpm["load"] = mpm["load-diff"]
        if load_zeroing:
            mpm["load"] = mpm["load"] - mpm["load"].values[0]
        mpm["stress"] = mpm["load"] / (0.06 - mpm["disp"])
    return mpm

sb1 = pd.read_csv("./data/sb1.csv")
sb2 = pd.read_csv("./data/sb2.csv")
sb3 = pd.read_csv("./data/sb3.csv")
sb4 = pd.read_csv("./data/sb4.csv")

plt.plot(sb1["x"].values,sb1["y"].values,label="SB1, SN = 401.5kPa")
plt.plot(sb2["x"].values,sb2["y"].values,label="SB2, SN = 289kPa")
plt.plot(sb3["x"].values,sb3["y"].values,label="SB3, SN = 184kPa")
plt.plot(sb4["x"].values,sb4["y"].values,label="SB4, SN = 72.5kPa")
def plot_mpm(folder,c):
    if os.path.isdir("../"+folder):
        mpm = get_load("../{}/disp.csv".format(folder))
        if len(mpm["load"]) > 0:
            l=plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=folder,marker=".",c=c)

colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
# plot_mpm("output-10.0d0_4.0_3_0.5_1.0_1000.0-401000.0","o")
plot_mpm("output-1.0d0_4.0_3_0.1_1.0_1000.0-401000.0",colours[0])
plot_mpm("output-1.0d0_4.0_3_0.1_1.0_1000.0-289000.0",colours[1])
plot_mpm("output-1.0d0_4.0_3_0.1_1.0_1000.0-184000.0",colours[2])
plot_mpm("output-1.0d0_4.0_3_0.1_1.0_1000.0-72000.0",colours[3])
#output-1.0d0_4.0_2_1.0_1.0_100.0-289000.0/
#output-1.0d0_4.0_2_1.0_1.0_1000.0-289000.0/
#output-10.0d0_4.0_2_1.0_1.0_100.0-289000.0/
#output-10.0d0_4.0_2_1.0_1.0_1000.0-289000.0/

# plot_mpm("../output-T_4.0_2_1.0_1.0_100.0-401000.0")
plt.xlabel("Displacement (mm)")
plt.ylabel("Load (N)")
plt.legend()
plt.show()
