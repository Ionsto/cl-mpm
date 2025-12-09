import matplotlib.pyplot as plt
import json
import re
import os
import pandas as pd

top_dir = "../"
output_regex = re.compile("output-*")
output_list = list(filter(output_regex.match,os.listdir(top_dir)))
output_list.sort()
if len(output_list) > 1:
    for i,out in enumerate(output_list):
        print("{}: {}".format(i,out))
    output_dir = "{}./{}/".format(top_dir,output_list[int(input())])
else:
    output_dir = "{}./{}/".format(top_dir,output_list[0])
df = pd.read_csv(output_dir+"timesteps.csv")
print(df)
energy_max = 1e-2
oobf_max = 1e-2
threshold = 0.1
with open(output_dir+"./settings.json") as f:
    json_settings = json.load(f)
    h = json_settings["RESOLUTION"]
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
    thresh_energy = [0,json_settings["DOMAIN-SIZE"][1]]
    thresh_energy = json_settings["CRITERIA-ENERGY"]
    thresh_oobf = json_settings["CRITERIA-OOBF"]
    thresh_hist = json_settings["CRITERIA-HIST"]

fig,ax = plt.subplots()
ax_damage = ax.twinx()
ax.set_yscale("log")
time = df["time"].values
ax_damage.plot(time,df["damage"].values,label="Damage",c="red",marker="x")
# plt.plot(df["time"].values,df["plastic"].values/df["plastic"].max() ,label="Plastic")
# plt.plot(df["time"].values,df["mass"].values/df["mass"].max() ,label="Mass")
# ax.plot(df["time"].values,df["energy"].values/df["work"].values,label="Energy",c="orange")
ax.plot(df["time"].values,df["oobf"].values,label="OOBF",c="green")
# ax.plot(df["time"].values,df["energy"].values,label="energy",c="orange")
quasi_point = None
for i in range(len(df)-1):
    if df["step-type"].iloc[i] == "DYNAMIC":
        if quasi_point == None:
            quasi_point = time[i]
    else:
        if quasi_point:
            plt.axvspan(quasi_point,time[i+1],alpha=0.25,color="black")
            quasi_point = None
if quasi_point:
    plt.axvspan(quasi_point,time[i+1],alpha=0.25,color="black")
# for i in range(len(df)-1):
#     if df["step-type"].iloc[i] != df["step-type"].iloc[i+1]:
#         ##Transition found
#         x = df["time"].values[i+1]
#         plt.axvline(x)
#         plt.text(x+1,0,'Transition to {}'.format(df["step-type"].values[i+1]),rotation=90)
# plt.axhline(threshold*(thresh_energy*thresh_hist)/thresh_energy,c="green",ls="--")
# plt.axhline(threshold*(thresh_energy/thresh_hist)/thresh_energy,c="red",ls="--")
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax_damage.get_legend_handles_labels()
ax_damage.legend(lines+lines2, labels+labels2, loc=0)
ax.set_xlabel("Time (s)")
ax.set_ylabel("Resiudal")
ax_damage.set_ylabel("Mass-damage")
plt.tight_layout()
plt.show()
