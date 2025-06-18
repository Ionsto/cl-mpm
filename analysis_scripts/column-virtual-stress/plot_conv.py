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
output_regex = re.compile("data-*")
output_list = list(filter(output_regex.match,os.listdir(top_dir)))
output_list.sort()
for i,out in enumerate(output_list):
    print("{}: {}".format(i,out))
data_h = []
data_e = []
for f in output_list:
    L = 50
    pressure = -1e4
    output_dir = "{}./{}".format(top_dir,f)
    x,y = f.split("-")[-1].split(".")[0].split("_")
    refine = int(x)
    mps = int(y)
    df = pd.read_csv(output_dir)
    y = df["Y"].values
    syy = df["SYY"].values
    syy_ref = df["SYY-REF"].values
    vp = df["VP"].values
    syy-syy_ref
    h = L / 2**refine
    data_h.append(h)
    data_e.append(np.linalg.norm(syy-syy_ref)*vp[0]/(np.sum(vp)))

data_h = np.array(data_h)
data_e = np.array(data_e)
plt.scatter(1/data_h,data_e,label="MPM")
plt.xlabel("1/h")
plt.ylabel("normalised error")
plt.xscale("log")
plt.yscale("log")
plt.show()
