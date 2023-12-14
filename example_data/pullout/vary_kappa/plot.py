import matplotlib.pyplot as plt
import pandas as pd
import os

files = os.listdir("csvs")

for i in files:
    mpm = pd.read_csv("csvs/{}".format(i))
    plt.plot(mpm["disp"],mpm["load"],label="Reaction - h/{}".format(i))
plt.xlabel("Displacement")
plt.ylabel("Load")
plt.legend()
plt.show()
