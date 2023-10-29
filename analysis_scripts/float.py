import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
plt.close("all")
for h in [50,25,10]:
    uri = "output/surface_position_{}.csv".format(h)
    if Path(uri).is_file():
        df = pd.read_csv(uri)
        plt.plot(df["Time (s)"],df["Surface position"]-300,label=h)
plt.axhline(10)
plt.legend()
plt.show()
