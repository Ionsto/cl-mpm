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



top_dir = "../../../"
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
# names = ["K constant","K updated","P elastic","P elastoplastic"]
names = ['output-K-0', 'output-K-UPDATED', 'output-P-ELASTIC', 'output-P-PLASTIC-SCALAR', 'output-P-PLASTIC-TANGENT']
for n,out in zip(names,output_list):
    output_dir = "{}./{}/".format(top_dir,out)
    df = pd.read_csv(output_dir+"conv.csv")
    c = None
    data_steps = []
    data_rates = []
    for name,group in df.groupby("step"):
        iters = group["iter"].values
        iters = iters - iters[0]
        oobf = group["oobf"].values
        data_steps.append(name)
        res = np.log10(oobf)
        # l = plt.plot(iters,res)
        m,b = np.polyfit(iters, res, 1)
        data_rates.append(m)
        # if c == None:
        #     l = plt.plot(iters,oobf,label=n)
        #     c = l[0].get_color()
        # else:
        #     plt.plot(iters,oobf,c=c)
    data_rates = np.array(data_rates)
    plt.plot(data_steps,-data_rates,label=n)
    # print(df["iters"].values[-1])

# thresh_scale = 1e-9
# plt.axhline(thresh_scale,c="green",ls="--")
# ax.set_ylim(bottom=0,top=thresh_scale_damage*2)
plt.xlabel("Load step")
plt.ylabel("Average convergence rate")
plt.yscale("log")
#plt.legend(["Aggregated","Non-aggregated"])
plt.legend()
plt.tight_layout()
plt.savefig("rate_comp.pdf")
plt.show()
