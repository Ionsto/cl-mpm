PDF_OUTPUT = False
import matplotlib as mpl
if PDF_OUTPUT:
    mpl.use('pdf')
else:
    mpl.use('Agg')
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

NO_OVERWRITE = False


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



import subprocess

plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
width = 3.487
height = width / 1.618

water_height = 236
xlim = [0,2600]
ylim = [0,500]
with open("output/settings.json") as f:
    json_settings = json.load(f)
    #print("Water level:{}".format(json_settings["OCEAN-HEIGHT"]))
    #print("Domain size:{}".format(json_settings["DOMAIN-SIZE"]))
    #water_height = json_settings["OCEAN-HEIGHT"]
    water_height = 0
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
    #ylim[0] = 20


ice_height = 200

plt.close("all")
output_dir = "./output/"
files = os.listdir(output_dir)
finalcsv = re.compile("sim.*\.vtk")
files_csvs = list(filter(finalcsv.match,files))
print("files: {}".format(files_csvs))
dt = 1e4/60
time = []
max_stress = []
damage = []
full_data = []

if not NO_OVERWRITE:
    subprocess.run("rm ./outframes/*", shell=True)

fig = plt.figure(figsize=(16,9),dpi=200)
def get_plot(i,fname):
    outname = "outframes/frame_{:05}.png".format(i)
    if NO_OVERWRITE and os.path.isfile(outname):
        return
    df = get_data(output_dir+"{}".format(fname))
    print("Plot frame {}".format(i),flush=True)
    ax = fig.add_subplot(111,aspect="equal")
    #df = full_data[i]
    patch_list=[]
    patch = Rectangle(xy=(0,0) ,width=xlim[1], height=water_height,color="blue")
    patch_sea = [patch]
    ps = PatchCollection(patch_sea)
    ax.add_collection(ps)

    for a_x, a_y,lx,ly,damage in zip(df["coord_x"],
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
    fig.colorbar(p,location="bottom",label="damage")

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.title("t = {:.2f}s".format(i * dt))
    plt.savefig("outframes/frame_{:05}.png".format(i))
    plt.clf()

def wrapper(x):
    get_plot(x[0],x[1])
if __name__ == '__main__':
    with Pool(8) as p:
        p.map(wrapper, enumerate(files_csvs))
cmd_str = "ffmpeg -y -framerate 60 -pattern_type glob -i 'outframes/*.png' -c:v libx264 -pix_fmt yuv420p out.mp4"
subprocess.run(cmd_str, shell=True)

