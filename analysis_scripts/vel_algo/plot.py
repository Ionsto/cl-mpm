import matplotlib as mpl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
from scipy.ndimage import uniform_filter1d
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

def get_data(filename):
    global data_name
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
    #damage = GetScalar("plastic_strain")
    #damage = GetScalar("damage")
    return pd.DataFrame({"coord_x":xy[:,0], "coord_y":xy[:,1],"lx":lx,"ly":ly,
                         "mass":GetScalar("mass"),
                         "vel_x":GetScalar("vel_x"),
                         "vel_y":GetScalar("vel_y")
                         })

full_data = []
for name in ["FLIP","PIC","BLEND"]:
    df = get_data("./data/sim_{}_5.0d-1.vtk".format(name))
    energy = 0.5 * df["mass"].values * (df["vel_x"].pow(2) + df["vel_y"].pow(2)).pow(0.5)
    full_data = np.hstack((full_data,energy))

bins = np.histogram(full_data, bins=2**12)[1] # Get the bin edges
for name in ["FLIP","PIC","BLEND"]:
    df = get_data("./data/sim_{}_5.0d-1.vtk".format(name))
    energy = 0.5 * df["mass"].values * (df["vel_x"].pow(2) + df["vel_y"].pow(2)).pow(0.5)
    plt.hist(energy,bins=bins,label=name,alpha=0.5)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()
