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



# top_dir = "./data/"
# output_regex = re.compile("data-*")
# output_list = list(filter(output_regex.match,os.listdir(top_dir)))
# output_list.sort()
# for i,out in enumerate(output_list):
#     print("{}: {}".format(i,out))
# output_dir = "{}./{}".format(top_dir,output_list[int(input())])


output_dir = "data.csv"
df = pd.read_csv(output_dir)
x = df["X"].values
ma = df["MA"].values
mii = df["M-LUMPED"].values
plt.ion()
plt.plot(x+0.5,ma,label="Aggregate")
plt.plot(x+0.5,mii,label="Lumped")
plt.xlabel("Displacement")
plt.ylabel("Condition number")
plt.yscale("log")
plt.xlim(0.5,1.5)
plt.gca().set_xticks(np.linspace(0.5,1.5,5))
plt.legend()
plt.show()
