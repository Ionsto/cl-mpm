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
    
water_height = 236
xlim = [0,2600]
ylim = [0,500]
with open("output_mici/settings.json") as f:
    json_settings = json.load(f)
    print("Water level:{}".format(json_settings["OCEAN-HEIGHT"]))
    print("Domain size:{}".format(json_settings["DOMAIN-SIZE"]))
    water_height = json_settings["OCEAN-HEIGHT"]
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
    ylim[0] = 20

for i in range(10,12):
    fname = "output_mici/sim_000{}.vtk".format(i)
    fig = plt.figure()
    df = get_data(fname)
    plot_data(df,xlim,ylim)
plt.show()
