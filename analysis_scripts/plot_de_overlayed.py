import matplotlib.pyplot as plt
import json
import pandas as pd
import re
import os

top_dir = "../"
output_regex = re.compile("output-*")
output_list = list(filter(output_regex.match,os.listdir(top_dir)))
output_list.sort()
for i,out in enumerate(output_list):
    print("{}: {}".format(i,out))


fig = plt.figure()
ax = fig.gca()
for outdir in output_list:
    output_dir = "{}./{}/".format(top_dir,outdir)
    df = pd.read_csv(output_dir+"timesteps.csv")
    time = df["time"].values
    step = df["steps"].values
    damage = df["damage"].values
    plt.plot(time,damage ,label=outdir,marker="x")

plt.xlabel("Time")
plt.ylabel("Mass-damage")

plt.legend()
plt.show()
