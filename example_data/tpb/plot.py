import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("load-disp.csv")
plt.plot(data["disp"],data["load"],label="Data")

scales = [1,1,1]
for i in ["0.5","1","2"]:
    mpm = pd.read_csv("conv/disp-{}.csv".format(i))
    plt.plot(-1e3*mpm["disp"],13e-3*mpm["load"],label="Reaction - h/{}".format(i))
for i in ["4"]:
    mpm = pd.read_csv("conv/disp-{}.csv".format(i))
    plt.plot(-1e3*mpm["disp"],13e-3*mpm["nload"],label="Reaction - h/{}".format(i))
plt.xlabel("Displacement (mm)")
plt.ylabel("Load (N)")
plt.legend()
plt.show()
