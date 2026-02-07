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
    L = 1
    output_dir = "{}./{}".format(top_dir,f)
    x,y = f.split("-")[-1][:-4].split("_")
    refine = float(x)
    mps = int(y)
    df = pd.read_csv(output_dir)
    # y = df["Y"].values
    # syy = df["SYY"].values
    # syy_ref = df["SYY-REF"].values
    # vp = df["VP"].values
    # syy-syy_ref
    h = L / 2**refine
    e = float(df["ERROR"].values[0].replace("d","e"))
    #e = np.linalg.norm(syy-syy_ref)*vp[0]/(np.sum(vp)))

    data_h.append(h)
    data_e.append(e)

data_h = np.array(data_h)
data_e = np.array(data_e)
plt.scatter(1/data_h,data_e,label="MPM")

m,b = np.polyfit(np.log(1/data_h),np.log(data_e),1)
print(m)

def plot_tri(offset,size):
    xoffset = 10**offset[0]
    yoffset = 10**offset[1]
    xsize = 1+size[0]
    ysize = 1+size[1]
    x = [xoffset*xsize,xoffset*xsize ,xoffset]
    y = [yoffset      ,yoffset/ysize ,yoffset]
    pos = np.vstack((x,y)).transpose()
    print(pos)
    t1 = plt.Polygon(pos,color="black",fill=False,lw=1)
    plt.gca().add_patch(t1)
    plt.text((xoffset*xsize)*1.1,yoffset/(ysize**0.6),size[1],size="x-small")
    # plt.plot(x,y)

#plot_tri([1,-3],[1,1])
plot_tri([0.5,-2.0],[1,1])

plt.xlabel("1/h")
plt.ylabel("normalised error")
plt.xscale("log")
plt.yscale("log")
#plt.tight_layout()
plt.show()
