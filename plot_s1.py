import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
plt.figure()
data_dir = "./output/"
files = os.listdir(data_dir)
dt = 1e-2
p = re.compile('simcsv.*\.csv') 
h5_files = [f for f in files if p.match(f)]
max_stress = []
damage = []
time = []
for f in h5_files:
    df = pd.read_csv(data_dir+f)
    max_stress.append(df["stress_1"].max())
    damage.append(df["damage"].max())
    time.append(dt * 100 * float(re.findall("\d+",f)[0]))

time = np.array(time)
max_stress = np.array(max_stress)
damage = np.array(damage)
df = pd.DataFrame({"time":time,"stress_1":max_stress,"damage":damage})
df.to_csv("notch_viscoelastic.csv")
plt.figure()
plt.plot(time,max_stress/1e6,label="stress 1")
plt.plot(time,damage,label="damage")
plt.legend()
plt.show()

