import matplotlib.pyplot as plt
import pandas as pd

#data = pd.read_csv("disp.csv")
#plt.plot(data["disp"],data["load"],label="Data")

for i in ["2","4","6"]:
    mpm = pd.read_csv("conv/disp-{}.csv".format(i))
    plt.plot(mpm["disp"],mpm["load"],label="Reaction - {}".format(i))
plt.legend()
plt.show()
