import pandas as pd
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt

from analysis_plot_kit.core import fig_generator

# read in data
# watertap\watertap\examples\flowsheets\property_comparison\RO\sal_outputs_results_RO_NaCl.csv
RO_NaCl = pd.read_csv(
    r"watertap\watertap\examples\flowsheets\property_comparison\RO\sal_outputs_results_RO_NaCl.csv"
)
RO_Sea = pd.read_csv(
    r"watertap\watertap\examples\flowsheets\property_comparison\RO\sal_outputs_results_RO_seawater.csv"
)
RO_Simple = pd.read_csv(
    r"watertap\watertap\examples\flowsheets\property_comparison\RO\sal_outputs_results_RO_simple.csv"
)

# find errors and create dataframe

results_col = ["LCOW", "SEC", "Membrane Area", "Operating Pressure"]
error_nacl = pd.DataFrame()
error_simple = pd.DataFrame()
for i in results_col:
    error_nacl[i] = pd.DataFrame(abs(RO_Sea[i] - RO_NaCl[i]) / RO_Sea[i]) * 100
    error_simple[i] = pd.DataFrame(abs(RO_Sea[i] - RO_Simple[i]) / RO_Sea[i]) * 100

# drop na from all dataframes
error_nacl = error_nacl.dropna().reset_index(drop=True)
error_simple = error_simple.dropna().reset_index(drop=True)

RO_Sea = RO_Sea.dropna()
RO_Sea = RO_Sea.reset_index(drop=True)
RO_NaCl = RO_NaCl.dropna()
RO_NaCl = RO_NaCl.reset_index(drop=True)
RO_Simple = RO_Simple.dropna()
RO_Simple = RO_Simple.reset_index(drop=True)

for i in results_col:
    figure = fig_generator.figureGenerator()
    figure.init_figure(width=6, height=2, ncols=2, sharex=True, sharey=True)
    vmin = min([min(error_nacl[i]), min(error_simple[i])])
    vmax = max([max(error_nacl[i]), max(error_simple[i])])
    figure.plot_map(
        xdata=RO_NaCl["Inlet Salinity"] * 100,
        ydata=RO_NaCl["# Water Recovery"],
        zdata=error_nacl[i],
        text=False,
        build_map=True,
        ax_idx=0,
        vmin=vmin,
        vmax=vmax,
    )
    figure.plot_map(
        xdata=RO_Simple["Inlet Salinity"] * 100,
        ydata=RO_Simple["# Water Recovery"],
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
        np.linspace(min(RO_NaCl["Inlet Salinity"]), max(RO_NaCl["Inlet Salinity"]), 5)
        * 100
    )
    yticks = np.linspace(
        min(RO_NaCl["# Water Recovery"]), max(RO_NaCl["# Water Recovery"]), 5
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
    figure.save_fig(name=i + "_RO")


for i in results_col:
    figure = fig_generator.figureGenerator()
    figure.init_figure(width=6, height=3, ncols=3, sharex=True, sharey=True)
    vmin = min([min(RO_NaCl[i]), min(RO_Simple[i]), min(RO_Sea[i])])
    vmax = max([max(RO_NaCl[i]), max(RO_Simple[i]), max(RO_Sea[i])])
    figure.plot_map(
        xdata=RO_Sea["Inlet Salinity"] * 100,
        ydata=RO_Sea["# Water Recovery"],
        zdata=RO_Sea[i],
        text=False,
        build_map=True,
        ax_idx=0,
        vmin=vmin,
        vmax=vmax,
    )
    figure.plot_map(
        xdata=RO_NaCl["Inlet Salinity"] * 100,
        ydata=RO_NaCl["# Water Recovery"],
        zdata=RO_NaCl[i],
        text=False,
        build_map=True,
        ax_idx=1,
        vmin=vmin,
        vmax=vmax,
    )
    figure.plot_map(
        xdata=RO_Simple["Inlet Salinity"] * 100,
        ydata=RO_Simple["# Water Recovery"],
        zdata=RO_Simple[i],
        text=False,
        build_map=True,
        ax_idx=2,
        vmin=vmin,
        vmax=vmax,
    )
    zticks = np.linspace(vmin, vmax, 6)
    figure.add_colorbar(zticks=zticks, zlabel=i)
    xticks = (
        np.linspace(min(RO_NaCl["Inlet Salinity"]), max(RO_NaCl["Inlet Salinity"]), 5)
        * 100
    )
    yticks = np.linspace(
        min(RO_NaCl["# Water Recovery"]), max(RO_NaCl["# Water Recovery"]), 5
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
    print(i + "_comparison_RO")
    figure.save_fig(name=i + "_comparison_RO")
