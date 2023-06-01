PDF_OUTPUT = False
import matplotlib as mpl
if PDF_OUTPUT:
    mpl.use('pdf')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
import os
import re


plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)

# width as measured in inkscape
width = 3.487
height = width / 1.618
ice_length = 2000

output_dir = "./output_notch/"

plt.figure(1)
#plt.figure(2)

output_dir = "./output_notch/"

files = os.listdir(output_dir)
finalcsv = re.compile("bottom_txx_.*.csv")
files_csvs = list(filter(finalcsv.match,files))
numbers = re.compile("\d*")
lengths = sorted([[int(y) for y in x if y][0] for x in map(numbers.findall,files_csvs)])
lengths_mpm = lengths
print("Notch lengths")
print(lengths)

max_stress = []
stress_pos = []
for dname in lengths:
    df = pd.read_csv(output_dir+"/bottom_txx_{}.csv".format(dname))
    max_stress.append(df["t_xx"])
    stress_pos.append(df["x"])
    plt.figure()
    plt.scatter(df["x"],df["t_xx"],label=dname)
plt.show()
width = 3.487
height = width / 1.618
