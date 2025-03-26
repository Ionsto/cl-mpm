PDF_OUTPUT = True
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

NO_OVERWRITE = True


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
    #damage = GetScalar("fric-contact")
    # damage = GetScalar("sig_xy")
    damage = GetScalar("q-undamaged")
    # damage = GetScalar("s_1") + GetScalar("pressure")
    #damage = GetScalar("plastic_strain")
    return pd.DataFrame({"coord_x":xy[:,0], "coord_y":xy[:,1]-offset,"lx":lx,"ly":ly,"damage":damage})

def get_data_all(folder,frame_number):
    regex = re.compile(r'sim_conv(_\d+)?_{}'.format(frame_number))
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

ice_height = 200

plt.close("all")
# output_dir = "/mnt/d/Temp/ham-float/"
for name,output_dir in zip(
        [
            "bench",
            "nobench",
            # "pressure",
            # "body",
            # "zero_pressure",
            # "zero_body"
         ],
    [
        "../../output-bench/",
        "../../output-nobench/",
        # "../../output-cryo-pressure/",
        # "../../output-cryo-body/",
        # "../../output-zero-pressure/",
        # "../../output-zero-body/"
     ]):
    if os.path.isdir(output_dir):
        print("file: {}".format(name))
        water_height = 236
        xlim = [0,2600]
        ylim = [0,500]
        with open(output_dir+"settings.json") as f:
            json_settings = json.load(f)
            #print("Water level:{}".format(json_settings["OCEAN-HEIGHT"]))
            #print("Domain size:{}".format(json_settings["DOMAIN-SIZE"]))
            #water_height = 0
            h = json_settings["RESOLUTION"]
            #offset = json_settings["OFFSET"]
            offset = 2*h
            water_height = json_settings["OCEAN-HEIGHT"]-offset
            xlim = [0,json_settings["DOMAIN-SIZE"][0]]
            ylim = [0,json_settings["DOMAIN-SIZE"][1]]
            aspect = 6
            ice_length = 400 * aspect
            xlim = [0,ice_length + 200]
            ylim = [0,450]
        files = os.listdir(output_dir)
        finalcsv = re.compile("sim_conv(_0+)?_\d+.vtk")
        files_csvs = list(filter(finalcsv.match,files))
        framenumber_regex = re.compile("\d+")
        files_csvs = list(map(lambda x: framenumber_regex.findall(x)[-1],files_csvs))
        files_csvs.sort(key=int)
        print("files: {}".format(files_csvs))
        dt = 1e4/60
        time = []
        max_stress = []
        damage = []
        full_data = []
        # if not NO_OVERWRITE:
        #     subprocess.run("rm ./outframes/*", shell=True)
        ratio = ylim[1]/xlim[1]
        print(ratio)
        #scale = 0.5
        scale = 1
        fig_x = scale*16
        fig_y = scale*16*ratio
        fig = plt.figure(figsize=(fig_x,fig_y),dpi=500)
        def get_plot(i,fname):
            df = get_data_all(output_dir,fname)
            print("Plot frame {}".format(i),flush=True)
            ax = fig.add_subplot(111,aspect="equal")
            patch_list=[]
            patch = Rectangle(xy=(ice_length,0) ,width=xlim[1], height=water_height,color="blue")
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

            cmap = cm.get_cmap(None, 10) 
            #cmap = cm.jet
            p = PatchCollection(patch_list, cmap=cmap, alpha=1)
            p.set_array(df["damage"])
            #p.set_clim([0,1.0])
            # p.set_clim([-1e5,1e5])
            ax.add_collection(p)
            # fig.colorbar(p,location="right",label="damage")

            plt.xlabel("(m)")
            plt.ylabel("(m)")
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            plt.margins(x=0,y=0)
            plt.tight_layout()
            plt.gcf().patch.set_alpha(0)
            cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
            plt.colorbar(p, cax=cax) # Similar to fig.colorbar(im, cax = cax)
            plt.savefig("poster_{}.pdf".format(name))
            #plt.savefig("poster_conv_{}.pdf".format(fname))
            # plt.show()
            plt.clf()
        def wrapper(x):
            get_plot(x,files_csvs[x])

        wrapper(len(files_csvs)-1)
