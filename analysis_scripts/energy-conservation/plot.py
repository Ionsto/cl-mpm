import matplotlib.pyplot as plt
import json
import pandas as pd
import re
import os

top_dir = "../../"
output_regex = re.compile("output-*")
output_list = list(filter(output_regex.match,os.listdir(top_dir)))

output_list.sort()
if len(output_list) > 1:
    for i,out in enumerate(output_list):
        print("{}: {}".format(i,out))
    output_dir = "{}./{}/".format(top_dir,output_list[int(input())])
else:
    output_dir = "{}./{}/".format(top_dir,output_list[0])


fig = plt.figure()
ax = fig.gca()
df = pd.read_csv(output_dir+"timesteps.csv")
time = df["time"].values
step = df["steps"].values
ke = df["ke"].values
se = df["se"].values
gpe = df["gpe"].values

# colour = l[0].get_color()
# plt.plot(time,se,c=colour)
# plt.plot(time,gpe,c=colour)
# plt.plot(time,ke + se - gpe,c=colour)
plt.plot(time,ke,label="kinetic energy")
plt.plot(time,se,label="strain energy")
plt.plot(time,gpe,label="gpe energy")
plt.plot(time,ke + se + gpe,label="energy delta")

plt.xlabel("Time")
plt.ylabel("Energy")

plt.legend()
plt.show()
