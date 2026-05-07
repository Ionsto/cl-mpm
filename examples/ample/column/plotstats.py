
PDF_OUTPUT = True
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1 import make_axes_locatable
# if PDF_OUTPUT:
#     mpl.use('pdf')
# else:
#     mpl.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
import matplotlib.ticker as plticker
import re
import os
import json
import numpy as np
import pandas as pd
import json
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
from multiprocessing import Pool
plt.style.use("seaborn-paper")
plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
plt.rcParams['figure.constrained_layout.use'] = True

width = 3.487
height = width / 1.618

fig = plt.figure(figsize=(width,height))

for f in ["stats.csv","stats_nr.csv"]:
    df = pd.read_csv(f)
    lstps = df["lstps"].values
    error = df["error"].values
    time = df["time"].values
    plt.plot(lstps,error,label=f)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.xlabel("Loadstep count")
plt.ylabel("Normalised stress error")

fig = plt.figure(figsize=(width,height))
for f in ["stats.csv","stats_nr.csv"]:
    df = pd.read_csv(f)
    lstps = df["lstps"].values
    time = df["time"].values
    plt.plot(lstps,time,label=f)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Loadstep count")
plt.ylabel("Solve time (s)")
plt.legend()
plt.show()

# plt.savefig("stats_convergence.pdf".format(output_name,i),dpi=1000,bbox_inches="tight",pad_inches=0)
