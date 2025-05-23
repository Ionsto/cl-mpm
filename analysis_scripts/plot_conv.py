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
for i,out in enumerate(output_list):
    print("{}: {}".format(i,out))
output_dir = "{}./{}/".format(top_dir,output_list[int(input())])


df = pd.read_csv(output_dir+"conv.csv")

sub_steps = 1
iters = df["iter"].values * sub_steps
step = df["step"].values
oobf = df["oobf"].values
energy = df["energy"].values
plastic = df["plastic"].values
plastic = plastic / np.max(plastic)
damage = df["damage"].values
damage = damage / np.max(damage)
fig = plt.figure()
ax = fig.gca()
#ax.plot(iters,oobf,label="OOBF")
#ax.plot(iters,energy,label="Energy")
ax.plot(iters,oobf,label="Residual")

thresh_scale = 1e-2
thresh_scale_damage = 1e-3
ax.axhline(thresh_scale,c="green",ls="--")
# ax.set_ylim(bottom=0,top=thresh_scale_damage*2)
ax.set_xlabel("Iterations")
ax.set_ylabel("Convergence criteria")
ax.set_yscale("log")
# ax.set_yscale("log")
ax_damage = ax.twinx()
ax_damage.plot(iters,damage,label="Damage",c="red")
ax_damage.plot(iters,plastic,label="Plastic",c="black")
ax_damage.set_ylim(bottom=0)
ax_damage.set_ylabel("Plastic strain evolution")

offset = 1
for i in range(len(df)-1):
    if step[i] != step[i+1]:
        ##Transition found
        x = iters[i+1]
        ax.axvline(x,color="black")
        ax.text(x+offset,-0.30,'Step {}'.format(step[i+1]),rotation=90)
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax_damage.get_legend_handles_labels()
ax_damage.legend(lines + lines2, labels + labels2, loc=0)
plt.tight_layout()
plt.show()
