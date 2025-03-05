import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv("throughputs.csv")
plt.plot(df["threads"],1/df["time"])
plt.show()
