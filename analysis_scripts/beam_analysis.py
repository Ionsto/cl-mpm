import matplotlib.pyplot as plt
import pandas as pd
df_noinc = pd.read_csv("output_rate_{}/final.csv".format("true"))
for dname in ["true","inc","logspin"]:
    df = pd.read_csv("output_rate_{}/final.csv".format(dname))
    plt.scatter(df["coord_x"],df["coord_y"],label=dname)
    dx = (df["coord_x"] - df_noinc["coord_x"]).pow(2)
    dy = (df["coord_y"] - df_noinc["coord_y"]).pow(2)
    e = (dx + dy).pow(0.5).sum()
    print("Error of {}: {}".format(dname,e))
plt.legend()
plt.show()
