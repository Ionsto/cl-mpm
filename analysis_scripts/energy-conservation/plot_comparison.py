import matplotlib.pyplot as plt
import json
import pandas as pd
import re
import os

top_dir = "../../"
output_regex = re.compile("output-*")
output_list = list(filter(output_regex.match,os.listdir(top_dir)))
output_list.sort()
# for i,out in enumerate(output_list):
#     print("{}: {}".format(i,out))


fig = plt.figure()
ax = fig.gca()
for outdir in output_list:
    output_dir = "{}./{}/".format(top_dir,outdir)
    df = pd.read_csv(output_dir+"timesteps.csv")
    time = df["time"].values
    step = df["steps"].values
    ke = df["ke"].values
    se = df["se"].values
    gpe = df["gpe"].values

    plt.plot(time,ke + se + gpe,label=outdir)
    # plt.plot(time,ke,label=outdir)

plt.xlabel("Time")
plt.ylabel("Energy")

plt.legend()
plt.show()
