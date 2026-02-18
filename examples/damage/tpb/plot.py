import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import os
import numpy as np
import re

from scipy import integrate

def extract_vals(f):
    output,refine,load = f.split("-")
    #refine = float(refine)
    return refine,float(load)

from scipy import integrate
def calculate_gf(disp,load):
    return integrate.trapz(load,disp)/(0.102*0.6*13e-3)

top_dir = "../../"#"./paper-1/damage-mc/"
#top_dir = "./paper-1/damage-mc/"
regex = re.compile(r'^output.*')
folders = list(filter(regex.search,os.listdir(top_dir)))

plt.figure(1)

# load_zeroing = True
load_zeroing = False
# load_combined = True
load_combined = False
load_clipping = False

data = pd.read_csv("load-disp.csv")

def get_load(filename):
    mpm = pd.read_csv(top_dir+filename)
    # data = pd.read_csv("load-disp.csv")
    # mpm = pd.read_csv("output/disp.csv")
    # plt.plot(data["disp"],data["load"],label="Data")
    # plt.plot(-1*mpm["disp"],0.012*mpm["load"],label="MPM")
    mpm["disp"] = mpm["disp"]
    mpm["load"] = mpm["load"]
    if load_clipping:
        mpm = mpm[mpm["disp"] >= 0.01e-3]
    if len(mpm["load"]) > 0:
        if load_zeroing:
            mpm["load"] = mpm["load"] - mpm["load"].values[0]
        mpm["stress"] = mpm["load"] / (0.06 - mpm["disp"])
    return mpm

print("GF experimental:",calculate_gf(1e-3*data["disp"],data["load"]))
folders.sort()
plt.plot(data["disp"].values,data["load"].values,label="Experimental data")

for i in folders:
    print("loading folder: ",i)
    mpm = get_load("./{}/disp.csv".format(i))
    if len(mpm["load"]) > 0:
        #if load_zeroing:
        #    mpm["load"] = mpm["load"] - mpm["load"].values[0]
        #l=plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
        plt.plot(-1e3*mpm["disp"].values,0.012*mpm["load"].values,label=i)
        # plt.plot(-1e3*mpm["disp"].values,0.012*mpm["load"].values,label=i)
        print("GF ",i," :",calculate_gf(-1*mpm["disp"],13e-3*mpm["load"]))
        # plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load-diff"].values,label=i,marker="x",c=l[0].get_color())
        #print("Shear modulus {}GPa".format(1e-9*mpm["load"].max()/mpm["disp"].values[mpm["load"].argmax()]))
        maxload = (1e-3/0.06)*mpm["load"].max()

plt.xlabel("Displacement (mm)")
plt.ylabel("Load (kN)")
plt.legend()
plt.show()
