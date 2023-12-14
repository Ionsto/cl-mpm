import matplotlib.pyplot as plt
import pandas as pd
import os

top_dir ="../../output_vary_kappa/" 

files = os.listdir(top_dir)

#files = ['output_1.0d0', 'output_0.5d0', 'output_0.75d0']

plt.figure()
kappa = []
gf = []
for i in files:
    mpm = pd.read_csv("{}{}/data.csv".format(top_dir,i))
    kappa.append(mpm["kappa"][0])
    gf.append(mpm["gf"][0])
plt.scatter(kappa,gf)
plt.legend()
plt.figure()
for i,name in zip(files,kappa):
    #plt.figure()
    #plt.title(i)
    mpm = pd.read_csv("{}{}/disp.csv".format(top_dir,i))
    plt.plot(mpm["disp"],mpm["load"],label="{}".format(name))
plt.xlabel("Displacement (mm)")
plt.ylabel("Load")
plt.legend()
plt.show()
