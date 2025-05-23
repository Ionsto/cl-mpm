import matplotlib.pyplot as plt
import json
import pandas as pd
output_dir = "../output/"
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

plt.plot(df["time"].values,df["damage"].values/df["damage"].max() ,label="Damage")
# plt.plot(df["time"].values,df["plastic"].values/df["plastic"].max() ,label="Plastic")
# plt.plot(df["time"].values,df["mass"].values/df["mass"].max() ,label="Mass")
plt.plot(df["time"].values,threshold*df["energy"].values / thresh_energy,label="Energy")
plt.plot(df["time"].values,threshold*df["oobf"].values   / thresh_oobf,label="OOBF")
for i in range(len(df)-1):
    if df["step-type"].iloc[i] != df["step-type"].iloc[i+1]:
        ##Transition found
        x = df["time"].values[i+1]
        plt.axvline(x)
        plt.text(x+1,0,'Transition to {}'.format(df["step-type"].values[i+1]),rotation=90)
plt.axhline(threshold*(thresh_energy*thresh_hist)/thresh_energy,c="green",ls="--")
plt.axhline(threshold*(thresh_energy/thresh_hist)/thresh_energy,c="red",ls="--")
plt.legend()
plt.show()
