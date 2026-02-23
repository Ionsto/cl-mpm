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
    i = np.argmax(load)
    print("Max at {}mm".format(disp[i]*1e3))
    return integrate.trapz(load[i:],disp[i:])

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
        plt.plot(mpm["disp"].values*1e3,mpm["load"].values,label=i)
        print("GF ",i," :",calculate_gf(mpm["disp"],mpm["load"]))

#d = 1
#h = 0.5
#rho = 2000
#g = 10
mu = 0.5
plt.axhline(mu,ls="--") 

plt.xlabel("Displacement (mm)")
plt.ylabel("Load (N)")
plt.legend()
plt.show()
