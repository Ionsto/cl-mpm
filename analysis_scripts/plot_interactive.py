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
    damage = GetScalar(data_name)
    #damage = GetScalar("plastic_strain")
    #damage = GetScalar("damage")
    return pd.DataFrame({"coord_x":xy[:,0], "coord_y":xy[:,1],"lx":lx,"ly":ly,"damage":damage})

def get_data_all(folder,frame_number):
    regex = re.compile(r'sim_conv(_0+)?_{}'.format(frame_number))
    files = list(filter(regex.search,os.listdir(folder)))
    subframes = [get_data(folder + "/" + f) for f in files]
    df = pd.concat(subframes)
    return df



import subprocess

plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
width = 3.487
height = width / 1.618


#output_regex = re.compile("output-\d+")
#output_list = list(filter(output_regex.match,os.listdir()))
#print(output_list)
# output_dir = "./{}/".format(output_list[int(input())])
output_dir = "./output/"
#output_dir = "./output/"
#xlim = [0.06,0.12+0.08]
#ylim = [0,0.10]
with open(output_dir+"settings.json") as f:
    json_settings = json.load(f)
    #print("Water level:{}".format(json_settings["OCEAN-HEIGHT"]))
    #print("Domain size:{}".format(json_settings["DOMAIN-SIZE"]))
    h = json_settings["RESOLUTION"]
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
    #ylim[0] = 20

water_height = 0

ice_height = 200

plt.close("all")

files = os.listdir(output_dir)
#Grab all files that are unique to first rank
finalcsv = re.compile("sim_conv(_0+)?_\d*\.vtk")
files_csvs = list(filter(finalcsv.match,files))

#finalcsv = re.compile("sim_0+_.*\.vtk")

framenumber_regex = re.compile("\d+")
files_csvs = list(map(lambda x: framenumber_regex.findall(x)[-1],files_csvs))
files_csvs.sort(key=int)
#files_csvs = list(map(lambda x: "sim_{}.vtk".format(x), files_csvs))
print("files: {}".format(files_csvs))
dt = 1e4/60
time = []
# timesteps = pd.read_csv(output_dir+"timesteps.csv")
max_stress = []
damage = []
full_data = []

def get_plot(i,fname):
    plt.clf()
    df = get_data_all(output_dir,fname)
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
    #p.set_clim([0,1.0])
    ax.add_collection(p)
    fig.colorbar(p,location="bottom",label=data_name)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    # plt.title("Frame {} - T {:.2f}s - {}".format(i,timesteps["time"].iloc[i],timesteps["step-type"].iloc[i]))
    # plt.savefig("outframes/frame_{:05}.png".format(i))
    # plt.clf()


current_frame = 0
max_frame = len(files_csvs)-1
def replot():
    get_plot(current_frame,files_csvs[current_frame])
    fig.canvas.draw()



data_name = "damage"
data_map = {
        "p":"plastic_strain",
        "d":"damage",
        "y":"damage-ybar-scaled",
        "u":"sig_xx",
        "o":"mpi-index",
        "x":"disp-x"
        }
print(data_map)
def on_press(event):
    global current_frame,data_name
    # print('press', event.key)
    sys.stdout.flush()
    # print(event.key)
    if event.key in data_map:
        data_name = data_map[event.key]
        replot()

    if event.key == 'i':
        data_name = "pressure"
        replot()
    if event.key == 's':
        current_frame = 0
        replot()
    if event.key == 'e':
        current_frame = max_frame
        replot()
    if event.key == 'right':
        current_frame = min(current_frame + 1,max_frame)
        replot()
    if event.key == 'left':
        current_frame = max(current_frame - 1,0)
        replot()
    if event.key == 'up':
        current_frame = min(current_frame + 10,max_frame)
        replot()
    if event.key == 'down':
        current_frame = max(current_frame - 10,0)
        replot()
    


print("Plotting {}, {}",current_frame,files_csvs[current_frame])

fig = plt.figure()#plt.figure(figsize=(16,9),dpi=200)
fig.canvas.mpl_connect('key_press_event', on_press)
replot()
plt.show()
#def wrapper(x):
#    get_plot(x[0],x[1])
#if __name__ == '__main__':
#    with Pool(8) as p:
#        p.map(wrapper, enumerate(files_csvs))
#cmd_str = "ffmpeg -y -framerate 60 -pattern_type glob -i 'outframes/*.png' -c:v libx264 -pix_fmt yuv420p out.mp4"
#subprocess.run(cmd_str, shell=True)

