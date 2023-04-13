import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
impor tos
import re


plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)

# width as measured in inkscape
width = 3.487
height = width / 1.618


files = os.listdir("output_notch")
finalcsv = re.compile("final.*.csv")
files_csvs = list(filter(finalcsv.match,files))
numbers = re.compile("\d*")
lengths = sorted([[int(y) for y in x if y][0] for x in map(numbers.findall,files_csvs)])
print("Notch lengths")
print(lengths)

#df_noinc = pd.read_csv("output_notch/final_{}.csv".format("true"))
#lengths = [10,20,30,40,50,80,100]
max_stress = []
stress_pos = []
for dname in lengths:
    df = pd.read_csv("output_notch/final_{}.csv".format(dname))
    plt.scatter(df["coord_x"],df["coord_y"],label=dname)
    s1 = df["eps"].max()
    max_stress.append(s1)
    stress_pos.append(df["coord_x"][df["eps"]==s1])
max_stress = np.array(max_stress)
stress_pos = np.array(stress_pos)

width = 3.487
height = width / 1.618

plt.legend()
plt.figure()
plt.plot(lengths,max_stress*1e-6,"-o")
plt.xlabel("Notch length ($m$)")
plt.ylabel("EPS ($MPa$)")

plt.gcf().subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
plt.gcf().set_size_inches(width, height)
plt.savefig("bench_stress.pdf")
plt.figure()
plt.plot(lengths,2000-stress_pos,"-o")
plt.xlabel("Notch length (m)")
plt.ylabel("S_1 distance from front")
plt.show()


