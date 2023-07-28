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


import subprocess

plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
width = 3.487
height = width / 1.618

ice_height = 400

water_height = 104.8

plt.close("all")
output_dir = "./output/"
files = os.listdir(output_dir)
finalcsv = re.compile("sim.*\.csv")
files_csvs = list(filter(finalcsv.match,files))
print("files: {}".format(files_csvs))
dt = 10e0
time = []
max_stress = []
damage = []
full_data = []
for i,dname in enumerate(files_csvs):
    df = pd.read_csv(output_dir+"{}".format(dname))
    s1 = df["stress_xx"].mean()
    max_stress.append(s1)
    damage.append(df["damage"].mean())
    time.append(i*dt)
    full_data.append(df)

subprocess.run("rm ./outframes/*", shell=True)

xlim = [0,5500]
ylim = [0,600]

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
    plt.title("t = {}s".format(i * dt))
    plt.savefig("outframes/frame_{:05}.png".format(i))
    #plt.gcf().subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
    #plt.gcf().set_size_inches(width, height)
    #plt.savefig("damage_states_{:05}.pdf".format(i))
    #p.set_clim([0,1])
    #plt.close("all")
    plt.clf()

#subprocess.run(["ls", "-l"]) 

# for complex commands, with many args, use string + `shell=True`:
#cmd_str = "ls -l /tmp | awk '{print $3,$9}' | grep root"
cmd_str = "ffmpeg -y -framerate 60 -pattern_type glob -i 'outframes/*.png' -c:v libx264 -pix_fmt yuv420p out.mp4"
subprocess.run(cmd_str, shell=True)
