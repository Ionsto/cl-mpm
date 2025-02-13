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
import subprocess

#output_dir = "../../output/"
output_dir = "/mnt/c/Temp/visco-ice/"
files = os.listdir(output_dir)
finalcsv = re.compile("sim(_0+)?_\d+.vtk")
files_csvs = list(filter(finalcsv.match,files))
framenumber_regex = re.compile("\d+")
files_csvs = list(map(lambda x: framenumber_regex.findall(x)[-1],files_csvs))
files_csvs.sort(key=int)
# files_csvs = list(map(lambda x: "sim_{}.vtk".format(x), files_csvs))
print("files: {}".format(files_csvs))

def get_data(filename,data_name):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    vtk_points = data.GetPoints()
    xyz3d = vtk_to_numpy( vtk_points.GetData() )
    xy = xyz3d[:,0:3]
    scalar_names = [reader.GetScalarsNameInFile(i) for i in range(0, reader.GetNumberOfScalarsInFile())]
    scalar_data = data.GetPointData()
    #scalar_names = scalar_data.GetArrayNames()
    def GetScalar(scalar_name):
        return vtk_to_numpy(scalar_data.GetArray(scalar_names.index(scalar_name)))
    lx = GetScalar("size_x")
    ly = GetScalar("size_y")
    lz = GetScalar("size_z")
    index = GetScalar("index")
    data = GetScalar(data_name)
    df = pd.DataFrame({"coord_x":xy[:,0],
                         "coord_y":xy[:,1],
                         "coord_z":xy[:,2],
                         "lx":lx,
                         "ly":ly,
                         "lz":lz,
                         "index":index,
                         "data":data})
    if index_apply_filter:
        df = df[df["index"]==index_filter]
    return df

def get_data_all(folder,frame_number,data_name):
    regex = re.compile(r'sim(_\d+)?_{}'.format(frame_number))
    files = list(filter(regex.search,os.listdir(folder)))
    # print(files)
    subframes = [get_data(folder + "/" + f,data_name) for f in files]
    # print("concat")
    df = pd.concat(subframes)
    return df

# os.makedirs("./csvs/")
def makefile(i,output_csv):
    f = files_csvs[i]
    print("File {} {}".format(i,f))
    df = get_data_all(output_dir,f,"viscosity")
    # print("start write")
    fname = output_csv+"./frame_{:05}.csv".format(i)
    with open(fname,"w") as fil:
        for i in range(len(df)):
            fil.write("{},{},{},{},{},{},{}\n".format(df["coord_x"].iloc[i],
                                                      df["coord_y"].iloc[i],
                                                      df["coord_z"].iloc[i],
                                                      df["lx"].iloc[i],
                                                      df["ly"].iloc[i],
                                                      df["lz"].iloc[i],
                                                      df["data"].iloc[i]))

index_apply_filter = True
index_filter = 0
makefile(90,"./csvs/")
index_filter = 1
makefile(90,"./csvs_1/")
# subprocess.run("rm ./csvs_ice/*", shell=True)
# if __name__ == '__main__':
#     with Pool(16) as p:
#         p.map(makefile, list(range(len(files_csvs))))
