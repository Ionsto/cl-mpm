import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("load-disp.csv")
plt.plot(data["disp"],data["load"],label="Data")

for i in [1,2]:
    mpm = pd.read_csv("conv/disp.csv")
    #plt.plot(-1e3*mpm["disp"],2*13e-3*mpm["load"],label="MPM")
    plt.plot(-1e3*mpm["disp"],2*0.5*13e-3*mpm["nload"],label="Reaction - {}".format(i))
plt.legend()
plt.show()
