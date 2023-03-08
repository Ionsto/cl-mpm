import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#df_noinc = pd.read_csv("output_notch/final_{}.csv".format("true"))
lengths = [10,20,30,40,50,80,100]
max_stress = []
stress_pos = []
for dname in lengths:
    df = pd.read_csv("output_notch/final_{}.csv".format(dname))
    plt.scatter(df["coord_x"],df["coord_y"],label=dname)
    s1 = df["stress_1"].max()
    max_stress.append(s1)
    stress_pos.append(df["coord_x"][df["stress_1"]==s1])
    # dx = (df["coord_x"] - df_noinc["coord_x"]).pow(2)
    # dy = (df["coord_y"] - df_noinc["coord_y"]).pow(2)
    # e = (dx + dy).pow(0.5).sum()
    # print("Error of {}: {}".format(dname,e))
plt.legend()
plt.figure()
plt.plot(lengths,max_stress,"-o")
plt.xlabel("Notch length (m)")
plt.ylabel("S_1")
plt.figure()
plt.plot(lengths,1000-np.array(stress_pos),"-o")
plt.xlabel("Notch length (m)")
plt.ylabel("S_1 distance from front")
plt.show()
