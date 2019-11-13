import numpy as np
import obspy
from obspy.taup import TauPyModel
from scipy import signal
import os
import matplotlib
from matplotlib import rc
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

# This script makes a contour plot of receiver function data in the time domain as a function of epicentral
# distance
# ------------------------------------------------------------------------------------------
# Last updated 11/12/2019 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------

# Read in receiver function data
data_directory = "/Users/aburky/IFILES/NETWORKS/"
ntwk = "IU"
stat = "PTCN"
gw = 1.0
# Undo time-shift if it existed
tshift = 10

# Construct path to receiver funcitons
rf_dir = data_directory + ntwk + "/" + stat + "/RFUNCS/GW" + ''.join(str(gw).split('.')) + "/"
rfs = obspy.read(rf_dir + '*.sac')
# Construct path to figure
fig_dir = data_directory + ntwk + "/" + stat + "/GRAPHICS/"
fig_name = ntwk + "." + stat + ".GW" + ''.join(str(gw).split('.')) + ".eps"

# Make sure figure directory exists
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Set matplotlib to use LaTeX fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Sort by epicentral distance
gcarcs = []
npts = []
delta = []
bad = []
k = 0
for i in range(0, len(rfs)):
    npts.append(i)
    npts[i] = rfs[i].stats.sac.npts
    delta.append(i)
    delta[i] = rfs[i].stats.sac.delta
    # Check if receiver function is anomalously long (make this automatic!)
    if rfs[i].stats.sac.e >= 150:
        bad.append(k)
        bad[k] = i
        k += 1
    # if rfs[i].stats.sac.user1 > 0.35 and rfs[i].stats.sac.e < 150:
    #     gcarcs.append(j)
    #     gcarcs[j] = [rfs[i].stats.sac.gcarc, i]
    #     j += 1

bad.sort(reverse=True)
for i in range(0, len(bad)):
    del npts[bad[i]], delta[bad[i]], rfs[bad[i]]

j = 0
for i in range(0, len(rfs)):
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
# Scale Factor
a = 10
for i in bins:
    binmin = i - 1
    binmax = i + 1
    idx = np.where(np.logical_and(sorted_gcarcs[:, 0] >= binmin, sorted_gcarcs[:, 0] < binmax))
    idxs = sorted_gcarcs[idx, 1].astype(int)
    rf_bin.append(j)
    rf_bin[j] = i + a*(np.sum(norm_rfs[idxs, :], axis=0)/len(idxs))
    if rf_bin[j].size == 0:
        rf_bin[j] = i + np.zeros((rf_bin[j].shape[0] + 1, rf_bin[j].shape[1]))
    j += 1

# Make moveout curves for predicted arrivals (using EQ depth of 100km?)
model = TauPyModel(model="iasp91")
phases = ["P", "P410s", "P660s"]
p410s = []
p660s = []
j = 0
ttbins = np.arange(25, 95, 1)
for i in ttbins:
    arrivals = model.get_travel_times(source_depth_in_km=100, distance_in_degree=i, phase_list=phases)
    p410s.append(j)
    p410s[j] = arrivals[1].time - arrivals[0].time
    p660s.append(j)
    p660s[j] = arrivals[2].time - arrivals[0].time
    j += 1

x = np.arange(-tshift, (npmax*dt)-tshift, dt)
# for i in range(0, len(shifted_rfs)):
#    plt.plot(x, shifted_rfs[i], 'k', linewidth=0.25)
# base = plt.gca().transData
# rot = transforms.Affine2D().rotate_deg(-90)
# trans = transforms.Affine2D().translate(0, 130)

rot = 0

# Make unrotated record section
if rot == 0:
    for i in range(0, len(rf_bin)):
        # for i, rf_b in enumerate(rf_bin[::-1]):
        plt.plot(x, rf_bin[i][0], 'k', linewidth=0.5)
        plt.fill_between(x, bins[i], rf_bin[i][0], where=rf_bin[i][0] > bins[i], facecolor='red')
    # Add travel time curves
    plt.plot(p410s, ttbins, 'k', linewidth=1)
    plt.plot(p660s, ttbins, 'k', linewidth=1)
    plt.ylim(29, 93)
    plt.xlim(0, (npmax*dt)-tshift)
    plt.xlabel('Time Relative to P Arrival (s)')
    plt.ylabel('Epicentral Distance $(^{\circ})$')
    plt.title("Receiver Function Data for {}.{}".format(ntwk, stat))
    plt.savefig(fig_dir + fig_name)
    plt.show()
else:
    for i in range(0, len(rf_bin)):
        plt.plot(rf_bin[i][0], x, 'k', linewidth=0.5)
        plt.fill_betweenx(x, bins[i], rf_bin[i][0], where=rf_bin[i][0] > bins[i], facecolor='red')
    plt.xlim(25, 95)
    plt.gca().invert_yaxis()
    plt.show()
