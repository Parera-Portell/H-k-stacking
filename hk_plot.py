# -*- coding: utf-8 -*-
"""
Joan A. Parera Portell

15/06/2019

This script plots H-K stacks.

Usage: python3 hk_plot.py [in file] [out plot] [Hmin] [Hmax] [kmin] [kmax]
 [bootstrap]

"""
import sys
import math as m
import numpy as np
import pandas as pd
import scipy.interpolate as sp
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as ep
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
rcParams['font.size'] = 15
rcParams['font.family'] = 'serif'

if len(sys.argv) < 8:
    print("Usage: python3 plthk.py [in file] [out plot] [Hmin] [Hmax] [kmin]\
 [kmax] [bootstrap]")
    sys.exit()

#--------------------------------Read input-----------------------------------#

infile = sys.argv[1]
outfile = sys.argv[2]
hmin = float(sys.argv[3])
hmax = float(sys.argv[4])
kmin = float(sys.argv[5])
kmax = float(sys.argv[6])
boot = sys.argv[7]

# Read files as tables and calculate covariance
suma = pd.read_csv(infile, sep=",", names=["h", "k", "s"])
bootstrap = pd.read_csv(boot, sep=",", header=0)
hste = float(bootstrap["Hste"])
kste = float(bootstrap["Kste"])
corr = float(bootstrap["CORR"])
covar = corr*hste*kste

#--------------------------Generate contour lines-----------------------------#
yi = np.linspace(kmin,kmax,50)
xi = np.linspace(hmin,hmax,50)
x,y = np.meshgrid(xi,yi)
zc = sp.griddata((suma["h"],suma["k"]),suma["s"],(x,y), method="linear")

#--------------------------Create 1 sigma ellipse-----------------------------#
# Locate data point of maximmum S
maxs_index = suma.loc[suma["s"]==suma["s"].max()].index[0]
# Get H and k (center of ellipse) at data point of maximum S
h = float(suma["h"][maxs_index]) 
k = float(suma["k"][maxs_index])
# Calculate covariance matrix
cov_mat = np.array([[hste**2, covar],[covar, kste**2]])
eigvalues, eigvectors = np.linalg.eig(cov_mat)
# Obtain axes from eigenvalues
a = 2*m.sqrt(eigvalues[0]) 
b = 2*m.sqrt(eigvalues[1])
# Calculation of ellipse tilt
tilt = m.degrees(0.5*m.atan(2*covar/((hste**2)-(kste**2))))

#--------------------------------Plotting-------------------------------------#
plt.figure(figsize=(5,4), tight_layout=True)
ax = plt.axes()
# Grid map
hkmap = ax.tripcolor(suma["h"],suma["k"],suma["s"], cmap="viridis", zorder=0)
ax.set_xlim((hmin,hmax))
ax.set_ylim((kmin,kmax))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
plt.colorbar(hkmap, aspect=25).set_label("Normalized stack amplitude",
            rotation=270, labelpad=20)
hkmap.set_clim(0,1)
# Contour map
ax.contour(x,y,zc, np.linspace(0.1,0.9,9), colors="k", linewidths=0.1, zorder=1)
ax.set_ylabel("$\kappa$", fontsize=16)
ax.set_xlabel("H (km)", fontsize=16)
# Error ellipse
ax.add_patch(ep((h,k), a, b, tilt, fill=False, linewidth=0.8, color="k", zorder=2))
# Central point
ax.scatter(h,k, 30, color="k", marker="*", zorder=3)
ax.annotate("H="+str(round(h,1))+" km, $\kappa$="+str(round(k,2)),
            xy=(0.95, 0.92), ha="right", xycoords='axes fraction', fontsize=14, 
            color="k")
filename = boot.split("/")[-1]
station = filename.split("_")[0]
ax.annotate(station, xy=(0.03, 0.03), ha="left", xycoords='axes fraction', 
            fontsize=14, color="k")
plt.savefig(outfile, bbox_inches="tight", dpi=300)
plt.show()
