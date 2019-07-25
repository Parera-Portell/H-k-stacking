# -*- coding: utf-8 -*-
"""
Joan A. Parera Portell

15/06/2019

This script plots H-K stacks.

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
    print("Usage: python3 hk_plot.py [in file] [out plot] [Hmin] [Hmax] [kmin]\
 [kmax] [bootstrap]")

#----------------------------------Input--------------------------------------#

infile = sys.argv[1]
outfile = sys.argv[2]
hmin = float(sys.argv[3])
hmax = float(sys.argv[4])
kmin = float(sys.argv[5])
kmax = float(sys.argv[6])
boot = sys.argv[7]

suma = pd.read_csv(infile, sep=",", names=["h", "k", "s"])
bootstrap = pd.read_csv(boot, sep=",", header=0)
hstd = float(bootstrap["Hstd"])
kstd = float(bootstrap["Kstd"])
corr = float(bootstrap["CORR"])
covar = corr*hstd*kstd

#--------------------------------Plotting-------------------------------------#
# Contour data
yi = np.linspace(kmin,kmax,150)
xi = np.linspace(hmin,hmax,150)
x,y = np.meshgrid(xi,yi)
zc = sp.griddata((suma["h"],suma["k"]),suma["s"],(x,y), method="linear")

# Error ellipse
maxs_index = suma.loc[suma["s"]==1].index[0] # find table index of maximum s
h = float(suma["h"][maxs_index]) # get H and k (center of ellipse)
k = float(suma["k"][maxs_index])
cov_mat = np.array([[hstd**2, covar],[covar, kstd**2]]) # covariance matrix
eigvalues, eigvectors = np.linalg.eig(cov_mat)
a = 2*m.sqrt(eigvalues[0]) 
b = 2*m.sqrt(eigvalues[1])
tilt = m.degrees(0.5*m.atan(2*covar/((hstd**2)-(kstd**2)))) # counterclockwise

# Plot h-k stack
plt.figure(figsize=(7,6))
ax = plt.axes()
hkmap = ax.tripcolor(suma["h"],suma["k"],suma["s"], cmap="viridis", zorder=0)
ax.set_xlim((hmin,hmax))
ax.set_ylim((kmin,kmax))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
plt.colorbar(hkmap, aspect=25, pad=0.04).set_label("Normalized stack amplitude",
            rotation=270, labelpad=20)
hkmap.set_clim(0,1)
ax.contour(x,y,zc, np.linspace(0.1,0.9,9), colors="k", linewidths=0.1, zorder=1)
ax.set_ylabel("$\kappa$", fontsize=25)
ax.set_xlabel("H (km)", fontsize=20)
ax.add_patch(ep((h,k), a, b, tilt, fill=False, linewidth=1.2, color="r", zorder=2))
ax.scatter(h,k, 70, color="k", marker="*", zorder=3)
ax.annotate("H="+str(round(h,1))+" km, $\kappa$="+str(round(k,2)),
            xy=(0.95, 0.92), ha="right", xycoords='axes fraction', fontsize=13, 
            color="w")
plt.savefig(outfile, bbox_inches="tight", dpi=300)
plt.show()
