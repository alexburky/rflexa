import numpy as np
import obspy
from scipy import signal
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
# import plotly.graph_objects as go
# from plotly.offline import iplot

# This script makes a contour plot of receiver function data in the time domain as a function of epicentral
# distance
# ------------------------------------------------------------------------------------------
# Last updated 11/11/2019 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------

# Read in receiver function data
rf_dir = "/Users/aburky/PycharmProjects/bermudaRFs/data/rfuncs/"
rfs = obspy.read(rf_dir + '*.sac')

# Sort by epicentral distance
gcarcs = []
npts = []
delta = []
j = 0
for i in range(0, len(rfs)):
    npts.append(i)
    npts[i] = rfs[i].stats.sac.npts
    delta.append(i)
    delta[i] = rfs[i].stats.sac.delta
    if rfs[i].stats.sac.user1 > 0.35:
        gcarcs.append(j)
        gcarcs[j] = [rfs[i].stats.sac.gcarc, i]
        j += 1

# print(np.sort(gcarcs, axis=0)[:, 0])

# Resample data for consistent plotting
npmax = np.max(npts)
dt = np.min(delta)
for i in range(0, len(rfs)):
    if rfs[i].stats.sac.npts != npmax:
        rfs[i].data = signal.resample(rfs[i].data, npmax)

x = np.arange(0, npmax*dt, dt)
y = np.sort(gcarcs, axis=0)[:, 0]
X, Y = np.meshgrid(x, y)
z = []
j = 0
for i in np.sort(gcarcs, axis=0)[:, 1]:
    z.append(j)
    z[j] = np.array(rfs[int(i)].data/np.max(rfs[int(i)].data))
    j += 1

# Plot
# fig = plt.subplot()
# CS = fig.contourf(X, Y, z)
plt.contourf(X, Y, z, 20, cmap='RdGy_r', levels=np.linspace(-0.05, 0.1, 20))
plt.xlim(30, 120)
plt.colorbar()
plt.show()
# cbar = fig.colorbar(CS)

# Using plotly (slow, but nice?)
# contour = go.Contour(x=x, y=y, z=z, contours=dict(start=0, end=0.1, size=0.01), line_smoothing=0.85)
# data = [contour]
# fig = go.Figure(data=data)
# fig.update_xaxes(range=[30, 100])
# fig.show(renderer='chrome')
