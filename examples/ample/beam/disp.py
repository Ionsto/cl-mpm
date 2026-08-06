import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


analytic_x = np.array(
[
   0
   ,0.083
   ,0.162
   ,0.235
   ,0.302
   ,0.494
   ,0.603
   ,0.670
   ,0.714
   ,0.744
   ,0.767
   ,0.785
   ,0.799
   ,0.811
])

analytic_y = np.array([
   0
   ,0.004
   ,0.016
   ,0.034
   ,0.056
   ,0.160
   ,0.255
   ,0.329
   ,0.388
   ,0.434
   ,0.472
   ,0.504
   ,0.531
   ,0.555
   ])

data = pd.read_csv("output-dr/disp.csv")
x_dr = data["x"].values[::-1]
y_dr = data["y"].values[::-1]
t_dr = data["type"].values[::-1]

data = pd.read_csv("output-nr/disp.csv")
x_nr = data["x"].values[::-1]
y_nr = data["y"].values[::-1]
t_nr = data["type"].values[::-1]

dpi = 200

width = (1920-(1080/2))/dpi
height =(1080)/dpi
fig = plt.figure(figsize=(width,height),dpi=dpi)
L = 10
for i in range(max(len(x_nr),len(x_dr))+1):
# for i in range(0,5):
    plt.plot(L*analytic_x,L*analytic_y,label="Analytic solution",color="black",ls="--")
    plt.plot(x_dr[:i],y_dr[:i],label="Dynamic Relaxation",lw=3)
    points_x = (x_dr[:i])[t_dr[:i]==1]
    points_y = (y_dr[:i])[t_dr[:i]==1]
    plt.scatter(points_x,points_y,color="C0",marker="x",s=100)

    plt.plot(x_nr[:i],y_nr[:i],label="Newton-Krylov",lw=3)
    points_x = (x_nr[:i])[t_nr[:i]==1]
    points_y = (y_nr[:i])[t_nr[:i]==1]
    plt.scatter(points_x,points_y,color="C1",marker="x",s=100)
    plt.xlim([-0.1,9])
    plt.ylim([-0.1,5.9])

    plt.xlabel("Horizontal tip displacement (m)")
    plt.ylabel("Vertical tip displacement (m)")
    plt.legend(loc="upper left")
    plt.savefig("disp_frames/disp_{:03d}.png".format(i))
    plt.clf()


