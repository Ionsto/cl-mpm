import pandas as pd
import matplotlib.pyplot as plt
for t in ["rate","jaumannrate","truesdallrate","true"]:
    df = pd.read_csv("deflection_{}.csv".format(t))
    plt.plot(df["x"],df["y"],label=t)
plt.legend()
plt.show()
