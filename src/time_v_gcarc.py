import numpy as np
import obspy
from scipy import signal
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

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

gcarcs = np.array(gcarcs)
sorted_gcarcs = gcarcs[gcarcs[:, 0].argsort()]

# Resample data for consistent plotting
npmax = np.max(npts)
dt = np.min(delta)
norm_rfs = []
for i in range(0, len(rfs)):
    if rfs[i].stats.sac.npts != npmax:
        rfs[i].data = signal.resample(rfs[i].data, npmax)
    norm_rfs.append(i)
    norm_rfs[i] = np.array(rfs[i].data/np.max(rfs[i].data))
norm_rfs = np.squeeze(norm_rfs)

# Construct vectors for plotting
shifted_rfs = []
for i in range(0, len(sorted_gcarcs)):
    shifted_rfs.append(i)
    shifted_rfs[i] = sorted_gcarcs[i][0] + \
        (rfs[int(sorted_gcarcs[i][1])].data/np.max(rfs[int(sorted_gcarcs[i][1])].data))
# Bin receiver functions in two-degree wide bins
bins = np.arange(30, 90, 1)
rf_bin = []
j = 0
for i in bins:
    binmin = i - 1
    binmax = i + 1
    idx = np.where(np.logical_and(sorted_gcarcs[:, 0] >= binmin, sorted_gcarcs[:, 0] < binmax))
    idxs = sorted_gcarcs[idx, 1].astype(int)
    rf_bin.append(j)
    rf_bin[j] = i + np.sum(norm_rfs[idxs, :], axis=0)/len(idxs)
    if rf_bin[j].size == 0:
        rf_bin[j] = i + np.zeros((rf_bin[j].shape[0] + 1, rf_bin[j].shape[1]))
    j += 1

x = np.arange(0, npmax*dt, dt)
# for i in range(0, len(shifted_rfs)):
#    plt.plot(x, shifted_rfs[i], 'k', linewidth=0.25)
for i in range(0, len(rf_bin)):
    plt.plot(x, rf_bin[i][0], 'k', linewidth=0.25)
plt.show()
