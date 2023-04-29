#import matplotlib as mpl
#mpl.use('pdf')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re


# plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
# plt.rc('xtick', labelsize=8)
# plt.rc('ytick', labelsize=8)
# plt.rc('axes', labelsize=8)

# width as measured in inkscape
width = 3.487
height = width / 1.618


files = os.listdir("output")
finalcsv = re.compile("sim.*.csv")
files_csvs = list(filter(finalcsv.match,files))
numbers = re.compile("\d*")
lengths = sorted([[int(y) for y in x if y][0] for x in map(numbers.findall,files_csvs)])
print("Notch lengths")
print(lengths)

#df_noinc = pd.read_csv("output_notch/final_{}.csv".format("true"))
#lengths = [10,20,30,40,50,80,100]
max_stress = []
stress_pos = []
files = []
for dname in files_csvs:
    df = pd.read_csv("output/{}".format(dname))
    #plt.scatter(df["coord_x"],df["coord_y"],label=dname)
    files.append(df)

df = files[0]
minx = df["coord_x"].min()
miny = df["coord_y"].min()
maxx = df["coord_x"].max()
maxy = df["coord_y"].max()




outside_ids = (df["coord_x"] == minx) | (df["coord_y"] == miny) | (df["coord_x"] == maxx) | (df["coord_y"] == maxy)
slicedf = df[outside_ids]
pos = [250,50]
angles = np.arctan2(slicedf["coord_y"]-pos[1],slicedf["coord_x"]-pos[0])
reorder = angles.argsort()

width = 3.487
height = width / 1.618

plt.figure()
plt.gca().set_aspect('equal')
#plt.scatter(slicedf["coord_x"].array[reorder],slicedf["coord_y"].array[reorder])
slices = [0,25,50,75,100]
for i in slices:
    if i < len(files):
        df = files[i]
        slicedf = df[outside_ids]
        x = slicedf["coord_x"].array[reorder]
        y = slicedf["coord_y"].array[reorder]
        plt.plot([*list(x),x[0]],[*list(y),y[0]])
        #plt.plot(slicedf["coord_x"].array[reorder][[-1,0]],slicedf["coord_y"].array[reorder][[-1,0]])
plt.show()
