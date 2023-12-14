import matplotlib.pyplot as plt
import pandas as pd
import os

files = os.listdir("csvs")

for i in files:
    mpm = pd.read_csv("csvs/{}".format(i))
    plt.plot(mpm["x"],mpm["d"],label="Reaction - h/{}".format(i))
plt.xlabel("Position (m)")
plt.ylabel("Damage (NA)")
plt.ylim([0,1])
plt.legend()
plt.show()
