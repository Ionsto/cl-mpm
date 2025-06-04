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



top_dir = "../"
output_regex = re.compile("output-*")
output_list = list(filter(output_regex.match,os.listdir(top_dir)))
output_list.sort()
fig = plt.figure()

# for i,out in enumerate(output_list):
#     print("{}: {}".format(i,out))

for out in output_list:
    output_dir = "{}./{}/".format(top_dir,out)
    df = pd.read_csv(output_dir+"conv.csv")
    df = df[df["step"]==1]
    x,y = out.split("-")[-1].split(".")[0].split("_")
    refine = int(x)
    mps = int(y)
    h = 1/ float((refine)**1)
    iters = df["iter"].values *h
    oobf = df["oobf"].values
    plt.plot(iters,oobf,label=out)

thresh_scale = 1e-9
plt.axhline(thresh_scale,c="green",ls="--")
# ax.set_ylim(bottom=0,top=thresh_scale_damage*2)
plt.xlabel("Iterations")
plt.ylabel("Convergence criteria")
plt.yscale("log")
plt.legend()
plt.show()
