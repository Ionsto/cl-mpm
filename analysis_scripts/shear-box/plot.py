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
    #*7.5e-3
    return 20*integrate.trapz(load,disp)/(0.06)

top_dir = "../../"#"./paper-1/damage-mc/"
#top_dir = "./paper-1/damage-mc/"
regex = re.compile(r'^output-.*')
folders = list(filter(regex.search,os.listdir(top_dir)))

plt.figure(1)

load_zeroing = True
# load_zeroing = False
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
    mpm["disp"] = -mpm["disp"]
    if load_clipping:
        mpm = mpm[mpm["disp"] >= 0.01e-3]
    if len(mpm["load"]) > 0:
        if load_zeroing:
            mpm["load"] = mpm["load"] - mpm["load"].values[0]
        mpm["stress"] = mpm["load"] / (0.06 - mpm["disp"])
    return mpm


data_name = []
data_gf = []
data_gf_target = []


folders.sort()
for i in folders:
    print("loading folder: ",i)
    mpm = get_load("./{}/disp.csv".format(i))
    if len(mpm["load"]) > 0:
        #if load_zeroing:
        #    mpm["load"] = mpm["load"] - mpm["load"].values[0]
        #l=plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
        plt.plot(-1e3*mpm["disp"].values,0.012*mpm["load"].values,label=i)
        raw_values = i.split("_")
        values = float(i.split("_")[3])
        print("GF ",i," :",calculate_gf(-1*mpm["disp"],13e-3*mpm["load"]))
        print("GF ratio ",i," :",calculate_gf(-1*mpm["disp"],13e-3*mpm["load"])/values)
        data_name.append(raw_values[0])
        data_gf.append(calculate_gf(-1*mpm["disp"],13e-3*mpm["load"]))
        data_gf_target.append(float(values))
        # plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load-diff"].values,label=i,marker="x",c=l[0].get_color())
        #print("Shear modulus {}GPa".format(1e-9*mpm["load"].max()/mpm["disp"].values[mpm["load"].argmax()]))
        maxload = (1e-3/0.06)*mpm["load"].max()
plt.xlabel("Displacement (mm)")
plt.ylabel("Load (kN)")
plt.legend()
gf_data = pd.DataFrame({"name":data_name,
                        "gf":data_gf,
                        "gf_target":data_gf_target,
                        }
                       )
#    np.array([data_name,data_gf,data_gf_target]).transpose(),columns=["name","gf","gf_target"])
plt.figure()
plt.figure()
for name,row in gf_data.groupby(by="name"):
    print(name)
    plt.scatter(row["gf_target"].values,row["gf"].values/row["gf_target"].values)
plt.show()
