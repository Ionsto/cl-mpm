import matplotlib as mpl
# if PDF_OUTPUT:
#mpl.use('pdf')
# else:
#mpl.use('Agg')
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
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
from multiprocessing import Pool
from matplotlib import pyplot, transforms

def get_data(filename):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    vtk_points = data.GetPoints()
    xyz3d = vtk_to_numpy( vtk_points.GetData() )
    xy = xyz3d[:,0:2]
    scalar_names = [reader.GetScalarsNameInFile(i) for i in range(0, reader.GetNumberOfScalarsInFile())]
    scalar_data = data.GetPointData()
    #scalar_names = scalar_data.GetArrayNames()
    def GetScalar(scalar_name):
        return vtk_to_numpy(scalar_data.GetArray(scalar_names.index(scalar_name)))
    lx = GetScalar("size_x")
    ly = GetScalar("size_y")
    damage = GetScalar("damage")
    return pd.DataFrame({"coord_x":xy[:,0], "coord_y":xy[:,1],"lx":lx,"ly":ly,"damage":damage})

def get_plot(i,fname):
    figscale = 0.3
    aspect = 1.0
                #16/9
    fig = plt.figure(figsize=(16*figscale,16*aspect*figscale),dpi=200)
    #outname = "outframes/frame_{:05}.png".format(i)
    #if NO_OVERWRITE and os.path.isfile(outname):
    #    return
    print(output_dir+"{}".format(fname))
    df = get_data(output_dir+"{}".format(fname))
    print("Plot frame {}".format(i),flush=True)
    ax = fig.add_subplot(111,aspect="equal")
    #df = full_data[i]
    patch_list=[]
    #patch = Rectangle(xy=(0,0) ,width=xlim[1], height=water_height,color="blue")
    #patch_sea = [patch]
    #ps = PatchCollection(patch_sea)
    #ax.add_collection(ps)

    for a_x, a_y,lx,ly,damage in zip(df["coord_x"]-(0 * 15.5),
                                     df["coord_y"],
                                     df["lx"],
                                     df["ly"],
                                     df["damage"]):
        patch = Rectangle(
            xy=(a_x-lx/2, a_y-ly/2) ,width=lx, height=ly)
        patch_list.append(patch)
    p = PatchCollection(patch_list, cmap=cm.jet, alpha=1)
    p.set_array(df["damage"])
    p.set_clim([0,1.0])
    ax.add_collection(p)
    #fig.colorbar(p,location="bottom",label="damage")

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticks([])
    ax.set_yticks([])
    #plt.title("t = {:.2f}s".format(i * dt))
    #plt.savefig("outframes/frame_{:05}.png".format(i))
    #plt.clf()

data_pre = pd.read_csv("data_pre.csv")
data_post = pd.read_csv("data_post.csv")

data_crack = pd.read_csv("cracks.csv")


plt.clf()
sync_number = 1
offset = 15.5
scale = 0.95

data_crack["x"] += 3
data_crack["y"] += 2.0

data_pre["x"] -= data_pre["x"].iloc[sync_number]
data_pre["y"] -= data_pre["y"].iloc[sync_number]
data_post["x"] -= data_post["x"].iloc[2] - 2
data_post["y"] -= data_post["y"].iloc[2]
sim_offset = [15.5-1.9 - 2,2]
scale = (15)/15.5
data_pre["y"] *=  scale
data_post["y"] *= scale
data_pre["y"] -= offset
data_post["y"] -=  offset


output_dir = "../../output/"
with open(output_dir+"settings.json") as f:
    json_settings = json.load(f)
    water_height = 0
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
xlim = [0,20]
    
files = os.listdir(output_dir)
finalcsv = re.compile("sim.*\.vtk")
files_csvs = list(filter(finalcsv.match,files))
print("files: {}".format(files_csvs))



rot = transforms.Affine2D().rotate_deg(5)
for i in [-1]:
    get_plot(i,files_csvs[i])
    #plt.plot(data_pre["x"] + sim_offset[0],-data_pre["y"] + sim_offset[1],label="Pre-failure geometry",linewidth=4.0,)
    #plt.plot(data_post["x"]+ sim_offset[0],-data_post["y"]+ sim_offset[1],label="Post-failure geometry",linewidth=4.0)
    plt.scatter(data_crack["x"],data_crack["y"],label="Tensile crack failure (Styles)",c="purple",marker="x")
    plt.legend()
    plt.savefig("plot.pdf")
plt.show()
