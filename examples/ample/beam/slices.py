import matplotlib.pyplot as plt
import mpmplotter
import mpmplotter.load
import mpmplotter.plot

plt.style.use("seaborn-paper")
plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
plt.rcParams['figure.constrained_layout.use'] = True

labels = {"nr":"Newton-Krylov","dr":"Dynamic Relaxation"}
for name in ["nr","dr"]:
    folder_name = "output-{}".format(name)
    data = mpmplotter.load.load_folder("./{}/".format(folder_name))
    offset = 2
    dpi = 100
    ratio = 1 #1.86 # 1.618
    width = 1
    height = width / ratio
    scale = 10
    fig = plt.figure(figsize=((1080/2)/dpi,(1080/2)/dpi),dpi=dpi)
    for i in range(len(data["frames"])):
        p = mpmplotter.plot.plot(data,i,colour_name="sig_xx")
        plt.xlim([0,12])
        plt.ylim([0,12])
        p.set_clim([-10e5,10e5])
        # fig.colorbar(p,location="right",label="damage")
        if name == "nr":
            ax = plt.gca()
            cax = ax.inset_axes([0.82, 0.05, 0.02, 0.5])
            cbar = fig.colorbar(p,cax=cax, orientation='vertical',label="Shear stress (Pa)")
        # cbar.ax.tick_params(labelsize=2)
        # cbar.set_label('Shear stress (MPa)', size=2)
        plt.title(labels[name])
        plt.xlabel("(m)")
        plt.ylabel("(m)")
        plt.savefig("outframes-{}/frame_{}.png".format(name,i))
        # plt.show()
        plt.clf()

        # 1080x
plt.close("all")
# plt.show()
