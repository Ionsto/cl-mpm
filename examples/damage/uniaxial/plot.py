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

top_dir = "../../../"
regex = re.compile(r'^output.*')
folders = list(filter(regex.search,os.listdir(top_dir)))
plt.figure(1)
def get_load(filename):
    mpm = pd.read_csv(top_dir+filename)
    mpm["disp"] = mpm["disp"]
    mpm["load"] = mpm["load"]
    return mpm


folders.sort()

for i in folders:
    print("loading folder: ",i)
    mpm = get_load("./{}/disp.csv".format(i))
    if len(mpm["load"]) > 0:
        plt.plot(mpm["disp"].values,mpm["load"].values,label=i)
        print("GF ",i," :",calculate_gf(mpm["disp"],mpm["load"]))
plt.xlabel("Displacement (mm)")
plt.ylabel("Load (kN)")
plt.legend()
plt.show()
