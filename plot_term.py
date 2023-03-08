import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.figure()
print("plotting")
df = pd.read_csv("output/terminus_position.csv")
plt.plot(df["Time (s)"],df["Terminus position"])
plt.legend()

plt.figure()
plt.plot(df["Time (s)"],np.gradient(df["Terminus position"],df["Time (s)"]))
#for f in [6,7]:
#    print("plotting {}".format(f));
#    df = pd.read_csv("output_1d{}/terminus_position.csv".format(f))
#    plt.plot(df["Time (s)"],df["Terminus position"]*(10**f),label=f)
#f = 6
#df = pd.read_csv("output_1d{}/terminus_position.csv".format(f))
#df["Terminus position"] -= df["Terminus position"][0]
#plt.plot(df["Time (s)"],df["Terminus position"],label=f)
#f = 7
#df = pd.read_csv("output_1d{}/terminus_position.csv".format(f))
#df["Terminus position"] -= df["Terminus position"][0]
#plt.plot(df["Time (s)"],df["Terminus position"]*10,label=f)
plt.show()

