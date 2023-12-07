import pandas as pd
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt

from analysis_plot_kit.core import fig_generator

# read in data
# watertap\watertap\examples\flowsheets\property_comparison\MVC\sal_outputs_results_MVC_NaCl.csv
MVC_NaCl = pd.read_csv(
    r"watertap\watertap\examples\flowsheets\property_comparison\MVC\sal_outputs_results_MVC_NaCl.csv"
)
MVC_Sea = pd.read_csv(
    r"watertap\watertap\examples\flowsheets\property_comparison\MVC\sal_outputs_results_MVC_seawater.csv"
)
MVC_Simple = pd.read_csv(
    r"watertap\watertap\examples\flowsheets\property_comparison\MVC\sal_outputs_results_MVC_simple.csv"
)

# find errors and create dataframe
results_col = ["LCOW", "SEC", "Evaporator area", "Compressor pressure ratio"]
error_nacl = pd.DataFrame()
error_simple = pd.DataFrame()
for i in results_col:
    error_nacl[i] = pd.DataFrame(abs(MVC_Sea[i] - MVC_NaCl[i]) / MVC_Sea[i]) * 100
    error_simple[i] = pd.DataFrame(abs(MVC_Sea[i] - MVC_Simple[i]) / MVC_Sea[i]) * 100

# drop na from all dataframes
error_nacl = error_nacl.dropna().reset_index(drop=True)
error_simple = error_simple.dropna().reset_index(drop=True)

MVC_Sea = MVC_Sea.dropna()
MVC_Sea = MVC_Sea.reset_index(drop=True)
MVC_NaCl = MVC_NaCl.dropna()
MVC_NaCl = MVC_NaCl.reset_index(drop=True)
MVC_Simple = MVC_Simple.dropna()
MVC_Simple = MVC_Simple.reset_index(drop=True)

for i in results_col:
    figure = fig_generator.figureGenerator()
    figure.init_figure(width=6, height=2, ncols=2, sharex=True, sharey=True)
    vmin = min([min(error_nacl[i]), min(error_simple[i])])
    vmax = max([max(error_nacl[i]), max(error_simple[i])])
    figure.plot_map(
        xdata=MVC_NaCl["Inlet Salinity"] * 100,
        ydata=MVC_NaCl["# Water Recovery"],
        zdata=error_nacl[i],
        text=False,
        build_map=True,
        ax_idx=0,
        vmin=vmin,
        vmax=vmax,
    )
    figure.plot_map(
        xdata=MVC_Simple["Inlet Salinity"] * 100,
        ydata=MVC_Simple["# Water Recovery"],
        zdata=error_simple[i],
        text=False,
        build_map=True,
        ax_idx=1,
        vmin=vmin,
        vmax=vmax,
    )
    zticks = np.linspace(vmin, vmax, 6)
    figure.add_colorbar(zticks=zticks, zlabel="% error, " + i)
    xticks = (
        np.linspace(min(MVC_NaCl["Inlet Salinity"]), max(MVC_NaCl["Inlet Salinity"]), 6)
        * 100
    )
    yticks = np.linspace(
        min(MVC_NaCl["# Water Recovery"]), max(MVC_NaCl["# Water Recovery"]), 7
    )
    xticks = np.around(xticks, decimals=1)
    yticks = np.around(yticks, decimals=2)
    figure.set_axis_ticklabels(
        xticklabels=xticks,
        yticklabels=yticks,
    )
    figure.set_fig_label(
        xlabel="Inlet Composition",
        ylabel="Recovery",
    )
    print(i)
    figure.save_fig(name=i)


for i in results_col:
    figure = fig_generator.figureGenerator()
    figure.init_figure(width=6, height=3, ncols=3, sharex=True, sharey=True)
    vmin = min([min(MVC_NaCl[i]), min(MVC_Simple[i]), min(MVC_Sea[i])])
    vmax = max([max(MVC_NaCl[i]), max(MVC_Simple[i]), max(MVC_Sea[i])])
    figure.plot_map(
        xdata=MVC_Sea["Inlet Salinity"] * 100,
        ydata=MVC_Sea["# Water Recovery"],
        zdata=MVC_Sea[i],
        text=False,
        build_map=True,
        ax_idx=0,
        vmin=vmin,
        vmax=vmax,
    )
    figure.plot_map(
        xdata=MVC_NaCl["Inlet Salinity"] * 100,
        ydata=MVC_NaCl["# Water Recovery"],
        zdata=MVC_NaCl[i],
        text=False,
        build_map=True,
        ax_idx=1,
        vmin=vmin,
        vmax=vmax,
    )
    figure.plot_map(
        xdata=MVC_Simple["Inlet Salinity"] * 100,
        ydata=MVC_Simple["# Water Recovery"],
        zdata=MVC_Simple[i],
        text=False,
        build_map=True,
        ax_idx=2,
        vmin=vmin,
        vmax=vmax,
    )
    zticks = np.linspace(vmin, vmax, 6)
    figure.add_colorbar(zticks=zticks, zlabel=i)
    xticks = (
        np.linspace(min(MVC_NaCl["Inlet Salinity"]), max(MVC_NaCl["Inlet Salinity"]), 6)
        * 100
    )
    yticks = np.linspace(
        min(MVC_NaCl["# Water Recovery"]), max(MVC_NaCl["# Water Recovery"]), 7
    )
    xticks = np.around(xticks, decimals=1)
    yticks = np.around(yticks, decimals=2)
    figure.set_axis_ticklabels(
        xticklabels=xticks,
        yticklabels=yticks,
    )
    figure.set_fig_label(
        xlabel="Inlet Composition",
        ylabel="Recovery",
    )
    print(i + "_comparison")
    figure.save_fig(name=i + "_comparison")
