PDF_OUTPUT = False
import matplotlib as mpl
if PDF_OUTPUT:
    mpl.use('pdf')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
import re
import os

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
    damage = GetScalar("damage_ybar")
    return pd.DataFrame({"coord_x":xy[:,0], "coord_y":xy[:,1],"lx":lx,"ly":ly,"damage":damage})

def plot_data(df,xlim,ylim):
    ax = fig.add_subplot(111,aspect="equal")
    patch_list=[]
    patch = Rectangle(xy=(0,0) ,width=xlim[1], height=water_height,color="blue")
    patch_sea = [patch]
    ps = PatchCollection(patch_sea)
    ax.add_collection(ps)

    for a_x, a_y,lx,ly in zip(df["coord_x"],
                                     df["coord_y"],
                                     df["lx"],
                                     df["ly"]):
        patch = Rectangle(xy=(a_x-lx/2, a_y-ly/2) ,width=lx, height=ly)
        patch_list.append(patch)
    p = PatchCollection(patch_list, cmap=cm.jet, alpha=1)
    p.set_array(df["damage"])
    #p.set_clim([0,1.0])
    p.set_clim([0,0.3e6])
    ax.add_collection(p)
    fig.colorbar(p,location="bottom",label="damage")

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


import subprocess

plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
width = 3.487
height = width / 1.618

ice_height = 400

water_height = 0.5 * 0.9 * ice_height

plt.close("all")
output_dir = "./output/"
files = os.listdir(output_dir)
finalcsv = re.compile("sim.*\.csv")
files_csvs = list(filter(finalcsv.match,files))
print("files: {}".format(files_csvs))
dt = 1e5
time = []
max_stress = []
damage = []
full_data = []
for i,dname in enumerate(files_csvs):
    #df = pd.read_csv(output_dir+"{}".format(dname))
    df = get_data(output_dir+"{}".format(dname))
    #s1 = df["stress_xx"].mean()
    #max_stress.append(s1)
    #damage.append(df["damage"].mean())
    time.append(i*dt)
    full_data.append(df)

subprocess.run("rm ./outframes/*", shell=True)

with open("output_mici/settings.json") as f:
    json_settings = json.load(f)
    print("Water level:{}".format(json_settings["OCEAN-HEIGHT"]))
    print("Domain size:{}".format(json_settings["DOMAIN-SIZE"]))
    water_height = json_settings["OCEAN-HEIGHT"]
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
    ylim[0] = 20


fig = plt.figure(figsize=(16,9),dpi=200)
for frame,i in enumerate(range(len(full_data))):
    print("Plot frame {}".format(i))
    ax = fig.add_subplot(111,aspect="equal")
    df = full_data[i]

    patch_list=[]
    patch = Rectangle(xy=(0,0) ,width=xlim[1], height=water_height,color="blue")
    patch_sea = [patch]
    ps = PatchCollection(patch_sea)
    ax.add_collection(ps)
    for a_x, a_y,lx,ly in zip(df["coord_x"],
                                     df["coord_y"],
                                     df["lx"],
                                     df["ly"]
                                     ):
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
    plt.title("t = {}s".format(i * dt))
    plt.savefig("outframes/frame_{:05}.png".format(i))
    plt.clf()

#subprocess.run(["ls", "-l"]) 

# for complex commands, with many args, use string + `shell=True`:
#cmd_str = "ls -l /tmp | awk '{print $3,$9}' | grep root"
cmd_str = "ffmpeg -y -framerate 30 -pattern_type glob -i 'outframes/*.png' -c:v libx264 -pix_fmt yuv420p out.mp4"
subprocess.run(cmd_str, shell=True)
