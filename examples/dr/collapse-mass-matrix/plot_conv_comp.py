import matplotlib as mpl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
import re
import os
import json
import numpy as np
import pandas as pd
import json
import sys
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
from multiprocessing import Pool



top_dir = "./data/"
output_regex = re.compile("output-*")
output_list = list(filter(output_regex.match,os.listdir(top_dir)))
output_list.sort()
fig = plt.figure()

# for i,out in enumerate(output_list):
#     print("{}: {}".format(i,out))

# for out in output_list:
#     output_dir = "{}./{}/".format(top_dir,out)
#     df = pd.read_csv(output_dir+"conv.csv")
#     iters = df["iters"].values
#     oobf = df["residual"].values
#     plt.plot(iters,oobf,label=out)

print(output_list)
names = ["K constant","K updated","P elastic","P elastoplastic"]
lss = ["-","-.",":","--"]
# lss = ["-","-","-","-"]
for n,out,ls in zip(names,output_list,lss):
    output_dir = "{}./{}/".format(top_dir,out)
    df = pd.read_csv(output_dir+"conv.csv")
    c = None
    for name,group in df.groupby("step"):
        iters = group["iters"].values
        oobf = group["residual"].values
        if c == None:
            l = plt.plot(iters,oobf,ls=ls,label=n)
            c = l[0].get_color()
        else:
            plt.plot(iters,oobf,ls=ls,c=c)
    print(out)
    print(df["iters"].values[-1])

ax = plt.gca()
x_1 = 11100
x_0 = 16820
y = 6e-10
ax.annotate('', xy=(x_0, y), xycoords='data',
xytext=(x_1,y), textcoords='data',
arrowprops=dict(arrowstyle="<->",
connectionstyle="bar", ec="k", shrinkA=5, shrinkB=5))
x_text = x_1 + 0.30 * (x_0 - x_1)
y_text = 5e-12
ax.annotate('lstp=10', xy=(x_text, y_text), xycoords='data')


x_1 = 9270
x_0 = 13740
y = 6e-10
ax.annotate('', xy=(x_0, y), xycoords='data',
xytext=(x_1,y), textcoords='data',
arrowprops=dict(arrowstyle="<->",
connectionstyle="bar", ec="k", shrinkA=5, shrinkB=5))
x_text = x_1 + 0.01 * (x_0 - x_1)
y_text = 8e-12
ax.annotate('lstp=9', xy=(x_text, y_text), xycoords='data')

# plt.axhline(thresh_scale,c="green",ls="--")
# ax.set_ylim(bottom=0,top=thresh_scale_damage*2)
plt.xlabel("Iterations")
plt.ylabel("Convergence criteria")
plt.yscale("log")
#plt.legend(["Aggregated","Non-aggregated"])
plt.legend()
plt.tight_layout()
plt.savefig("conv_comp.pdf")
plt.show()
